#!/bin/bash
#SBATCH -J 08_LD_diagnosis
#SBATCH -A naiss2025-5-732
#SBATCH -p shared
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-user=253721079@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/08_LD_diagnosis_%j.out
#SBATCH -e logs/08_LD_diagnosis_%j.err

set -euo pipefail

WORK=/cfs/klemming/projects/snic/snic2020-6-176/chen2017112/jingz-transfer/GS-904/GS-reviewer

D="${WORK}/08_imputed_GWAS/05_chip_imputed_LD_r2/outer_01/inner_01/chr_01"

LD="${D}/chip_imputed_r2.ld.gz"

[[ -s "${LD}" ]] || {
    echo "ERROR: Missing LD file" >&2
    exit 1
}

gzip -t "${LD}"

gzip -dc "${LD}" | awk \
-v chipfile="${D}/chip.snps" \
-v impfile="${D}/imputed.snps" '
BEGIN {
    while ((getline s < chipfile) > 0) {
        chip[s] = 1
    }
    close(chipfile)

    while ((getline s < impfile) > 0) {
        imp[s] = 1
        nimp++
    }
    close(impfile)
}

NR == 1 {
    for (j = 1; j <= NF; j++) {
        col[$j] = j
    }

    if (!("SNP_A" in col) ||
        !("SNP_B" in col) ||
        !("BP_A" in col) ||
        !("BP_B" in col) ||
        !("R2" in col)) {
        print "ERROR: Unexpected LD header" > "/dev/stderr"
        bad = 1
        exit 1
    }
    next
}

{
    a = $(col["SNP_A"])
    b = $(col["SNP_B"])
    r = $(col["R2"]) + 0

    snp = ""
    ref = ""

    if ((a in chip) && (b in imp)) {
        snp = b
        ref = a
    }
    else if ((b in chip) && (a in imp)) {
        snp = a
        ref = b
    }

    if (snp != "") {
        crosspairs++

        if (!(snp in maxr) || r > maxr[snp]) {
            maxr[snp] = r
            bestchip[snp] = ref

            distance = $(col["BP_A"]) - $(col["BP_B"])
            if (distance < 0) distance = -distance

            bestdistance[snp] = distance
        }
    }
}

END {
    if (bad || NR < 2) exit 1

    for (snp in imp) {
        if (!(snp in maxr)) {
            no_ld++
            continue
        }

        nmatched++
        r = maxr[snp]
        sum += r

        if (r < 0.2) n01_02++
        else if (r < 0.5) n02_05++
        else if (r < 0.8) n05_08++
        else if (r < 0.95) n08_095++
        else n095_1++

        if (r >= 0.2) above02++
        if (r >= 0.5) above05++
        if (r >= 0.8) above08++

        if (examples < 10) {
            example[++examples] = snp
        }
    }

    print "================================"
    print "Maximum chip-imputed LD analysis"
    print "================================"
    print "Total imputed SNPs:", nimp
    print "Imputed SNPs with r2 >= 0.1:", nmatched
    print "Imputed SNPs without r2 >= 0.1:", no_ld+0
    print "Chip-imputed pairs:", crosspairs+0

    if (nmatched > 0) {
        printf "Mean maximum r2: %.6f\n", sum/nmatched
    }

    print ""
    print "Maximum r2 distribution:"
    print "0.10 <= r2 < 0.20:", n01_02+0
    print "0.20 <= r2 < 0.50:", n02_05+0
    print "0.50 <= r2 < 0.80:", n05_08+0
    print "0.80 <= r2 < 0.95:", n08_095+0
    print "0.95 <= r2 <= 1.00:", n095_1+0

    print ""
    print "Imputed SNPs with max r2 >= 0.2:", above02+0
    print "Imputed SNPs with max r2 >= 0.5:", above05+0
    print "Imputed SNPs with max r2 >= 0.8:", above08+0

    print ""
    print "Examples:"
    print "Imputed_SNP Chip_SNP Max_R2 Distance_bp"

    for (j = 1; j <= examples; j++) {
        s = example[j]
        print s, bestchip[s], maxr[s], bestdistance[s]
    }
}
'
