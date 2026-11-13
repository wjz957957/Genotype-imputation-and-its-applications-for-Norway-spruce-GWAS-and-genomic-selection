# filter indel
# vcftools --gzvcf $1   --remove-indels --recode --recode-INFO-all --stdout | gzip -c > $fout".indel.recode.vcf.gz"
# filter biallelic
# vcftools --gzvcf $fout".indel.recode.vcf.gz"  --recode-INFO-all --max-alleles 2 --min-alleles 2 --recode  --stdout | gzip -c > $fout".indel.bia.recode.vcf.gz"

# 1)  only keep biallic and remodel indel, in addition, to remove SNPs that are within 10 bp of indels, we add the flag --SnpGap 10
# bcftools filter -e 'AC==0 || AC==AN'  --SnpGap 5  $1 | bcftools view -m2 -M2 -v snps   -s ^P_sitchensis,P_engelmannii,P_glauca --threads 4   -O z -o $fout".indel.bia.recode.vcf.gz"

# 2) DP 15000 is based on the first 100Mb, chr01_1, 3sd
# bcftools filter  -e 'QUAL<20 || DP>15000'   $fout".indel.bia.recode.vcf.gz"  --threads 4   -O z -o $fout".indel.bia.Qual20.DP.recode.vcf.gz"


# filter Missing 0.3, MAF
# bcftools filter -e 'F_MISSING > 0.3 || MAF <=0.0025 ' -O z -o $fout".indel.bia.miss.recode.vcf.gz"  $fout".indel.bia.recode.vcf.gz"


# rm $fout".indel.bia.recode.vcf.gz"

# get tabuate postion information,  space between \t and % is not necessary? yes
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $fout".indel.bia.Qual20.DP.recode.vcf.gz"  > pos.tabulated.txt

# rm $fout".indel.bia.Qual20.DP.recode.vcf.gz"

# awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$4}' pos.tabulated.txt > out.pos.tabulated.txt
# rm pos.tabulated.txt
# cp out.pos.tabulated.txt pos.tabulated.txt
# rm out.pos.tabulated.txt