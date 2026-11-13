#!/usr/bin/env python3
"""Choose QTN counts by SIGNED Pearson r with a two-stage fallback rule.

Priority 1: consider only nonzero QTN counts with r_k > r_0 in EACH of the
three inner folds, then choose the count with highest mean fold signed r.
Priority 2: when NO nonzero count passes, choose the count with the highest
mean fold signed r from ALL 20 counts, INCLUDING extra_0. If the baseline
ties for highest mean signed r, choose extra_0.

Never apply abs(r) or pooled Pearson for model selection. All comparisons
use the paired extra_0 baseline in the SAME inner fold. No outer-test is read.
"""
import argparse
import csv
import math
import os
from pathlib import Path
from collections import defaultdict

COUNTS=(0,1,2,3,4,5,10,20,30,40,50,75,100,150,200,300,400,500,600,700)


def check(cond,msg):
    if not cond:
        raise ValueError(msg)


def read_tsv(path):
    p=Path(path)
    check(p.is_file() and p.stat().st_size>0,f'missing {p}')
    with p.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def publish(path, rows, header):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp')
    with tmp.open('w', newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=header,delimiter='\t',extrasaction='ignore')
        wr.writeheader(); wr.writerows(rows)
    os.replace(tmp,p)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--work',required=True)
    ar=ap.parse_args()
    root=Path(ar.work)/'Cross_family/results/08_inner_qtn_rrblup'
    out=Path(ar.work)/'Cross_family/results/08_inner_qtn_selection'
    out.mkdir(parents=True,exist_ok=True)
    flag=out/'complete.flag'
    if flag.exists(): flag.unlink()
    rows_by_group={}; fold_rows=[]; n_preds=0
    all_pred_id_by_fold={}
    for o in range(1,6):
        for i in range(1,4):
            for t in ('BB','H6'):
                key=(f'{o:02d}',f'{i:02d}',t)
                p=root/f'outer_{o:02d}'/f'inner_{i:02d}'/t
                f=p/'complete.flag'
                check(f.is_file() and 'status=PASS' in f.read_text().splitlines() and
                      'baseline_reproduced=PASS' in f.read_text().splitlines(),f'flag invalid {f}')
                m=read_tsv(p/'metrics.tsv')
                check(len(m)==20, f'expected 20 metric rows: {p}')
                by_k={}
                nvalid=None
                for row in m:
                    check((row['outer'],row['inner'],row['trait'])==key,
                          f'metric fold/trait mismatch {p}')
                    k=int(row['n_qtn']); check(k in COUNTS and k not in by_k,f'duplicate/unexpected QTN count {p}')
                    check(row['status']=='PASS',f'model FAIL {p} extra_{k}')
                    for fld in ('pearson_r','RMSE','MAE','bias','delta_r','delta_RMSE'):
                        check(math.isfinite(float(row[fld])), f'nonfinite {fld}: {p} extra_{k}')
                    nv=int(row['n_valid']);check(nv>=50,'too few validation samples')
                    check(nvalid is None or nv==nvalid, f'inconsistent valid sample size {p}')
                    nvalid=nv;by_k[k]=row
                check(set(by_k)==set(COUNTS),f'incomplete QTN grid in {p}')
                base=float(by_k[0]['pearson_r'])
                check(abs(float(by_k[0]['delta_r'])) < 1e-12,'baseline delta_r nonzero')
                for k in COUNTS:
                    check(abs(float(by_k[k]['delta_r'])-(float(by_k[k]['pearson_r'])-base))<1e-8,
                          f'incorrect paired Pearson delta {p} extra_{k}')
                # Audit predictions for every candidate and model (distinct samples).
                pred=p/'validation_predictions.tsv'
                check(pred.is_file() and pred.stat().st_size>0, f'missing predictions {pred}')
                ids=defaultdict(set);observed={}
                with pred.open(newline='') as fh:
                    for row in csv.DictReader(fh,delimiter='\t'):
                        k=int(row['n_qtn']); sid=row['sample_id']
                        check(k in COUNTS and sid not in ids[k], f'duplicate prediction {pred} {k} {sid}')
                        check((row['outer'],row['inner'],row['trait'])==key,f'prediction mismatched {pred}')
                        yy=float(row['observed']); pr=float(row['predicted'])
                        check(math.isfinite(yy) and math.isfinite(pr),f'nonfinite prediction {pred}')
                        check(sid not in observed or abs(observed[sid]-yy)<1e-9,f'inconsistent observed {pred}')
                        observed[sid]=yy;ids[k].add(sid);n_preds+=1
                check(all(len(ids[k])==nvalid for k in COUNTS),f'incomplete prediction model coverage {pred}')
                check(all(ids[k]==ids[0] for k in COUNTS),f'differing valid sample subsets among models {pred}')
                all_pred_id_by_fold[key]=ids[0]
                rows_by_group[key]=by_k
                fold_rows.extend(m)
    ranking=[]; chosen=[]
    for o in range(1,6):
        for t in ('BB','H6'):
            outer=f'{o:02d}'
            observed_valid_sets=[all_pred_id_by_fold[(outer,f'{i:02d}',t)] for i in range(1,4)]
            concatenated=set()
            for z in observed_valid_sets:
                check(not concatenated.intersection(z),f'duplicate OOF validation sample: outer_{outer}/{t}')
                concatenated.update(z)
            outer_train_path=Path(ar.work)/'Cross_family/cv_samples_v1'/f'outer_{outer}'/'outer_train.samples'
            check(outer_train_path.is_file(),f'missing outer-training sample list: {outer_train_path}')
            full_train=[line.strip() for line in outer_train_path.read_text().splitlines()]
            check(len(full_train)==len(set(full_train)) and concatenated==set(full_train),
                  f'OOF validation coverage mismatch in outer_{outer}/{t}')
            folds=[rows_by_group[(outer,f'{i:02d}',t)] for i in range(1,4)]
            candidate=[]
            base_values=[float(fold[0]['pearson_r']) for fold in folds]
            for k in COUNTS:
                scores=[float(fold[k]['pearson_r']) for fold in folds]
                rmses=[float(fold[k]['RMSE']) for fold in folds]
                deltas=[s-b for s,b in zip(scores,base_values)]
                positive=sum(d>0 for d in deltas)
                eligible=k>0 and positive==3
                candidate.append(dict(outer=outer,trait=t,n_qtn=k,
                  mean_fold_pearson_r=sum(scores)/3,mean_fold_RMSE=sum(rmses)/3,
                  mean_delta_r=sum(deltas)/3,min_delta_r=min(deltas),
                  inner_01_delta_r=deltas[0],inner_02_delta_r=deltas[1],inner_03_delta_r=deltas[2],
                  n_folds_positive=positive,eligible='YES' if eligible else 'NO'))
            eligible_rows=[r for r in candidate if r['eligible']=='YES']
            mode=('THREE_OF_THREE_SIGNED_IMPROVEMENT' if eligible_rows else
                  'FALLBACK_MAX_MEAN_SIGNED_R')
            pool=eligible_rows if eligible_rows else candidate
            # Primary criterion = mean signed Pearson r, not mean abs(r).
            # On fallback, baseline wins an EXACT tie at the maximum mean r.
            # Nonbaseline ties use mean RMSE, then the smaller QTN count.
            best=sorted(pool,key=lambda r:(-r['mean_fold_pearson_r'],
                                           0 if r['n_qtn']==0 else 1,
                                           r['mean_fold_RMSE'],r['n_qtn']))[0]
            for r in candidate:
                r['selected']='YES' if r['n_qtn']==best['n_qtn'] else 'NO'
                r['selection_mode']=mode
                r['selection_rule']='all3_signed_improvement_else_max_mean_signed_Pearson'
                r['n_eligible_positive_counts']=len(eligible_rows)
            ranking.extend(candidate)
            chosen.append(best)
            print(f"outer_{outer}/{t}: selected extra_{best['n_qtn']} "
                  f"mean fold SIGNED r={best['mean_fold_pearson_r']:.6f} "
                  f"3/3 positive candidates={len(eligible_rows)} mode={mode}")
    ranking.sort(key=lambda r:(r['outer'],r['trait'],r['n_qtn']))
    selected_cols=['outer','trait','n_qtn','mean_fold_pearson_r','mean_fold_RMSE',
                   'mean_delta_r','min_delta_r',
                   'inner_01_delta_r','inner_02_delta_r','inner_03_delta_r',
                   'n_folds_positive','eligible','selection_mode',
                   'selection_rule','n_eligible_positive_counts']
    ranked_cols=list(ranking[0].keys())
    publish(out/'all_qtn_inner_fold_metrics.tsv',fold_rows,list(fold_rows[0]))
    publish(out/'qtn_count_comparison.tsv',ranking,ranked_cols)
    publish(out/'selected_qtn_counts.tsv',chosen,selected_cols)
    flag_tmp=flag.with_suffix('.flag.tmp')
    flag_tmp.write_text('version=3\nfold_traits=30\nmodel_metric_rows=600\nselection_groups=10\n'
                        f'prediction_rows={n_preds}\n'
                        'selection_rule=ALL_THREE_SIGNED_POSITIVE_ELSE_MAX_MEAN_SIGNED_R\n'
                        'ranking=MEAN_FOLD_SIGNED_PEARSON\n'
                        'absolute_value_selection=NO\n'
                        'outer_test_used=NO\nstatus=PASS\n')
    os.replace(flag_tmp,flag)
    print(f'OVERALL PASS: 30 fold traits, 600 models, 10 QTN decisions, {n_preds} predictions')

if __name__=='__main__':
    main()
