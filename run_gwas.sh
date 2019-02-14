#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q HighMemLongterm.q,LowMemLongterm.q,LowMemShortterm.q
#$ -t 1:20

module add bioinformatics/plink2/1.90b3.38
module add bioinformatics/R/3.3.3
module add compilers/gcc/6.2.0

# qsub -v file=sim5,bfile=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/sample_overlap/raw/ukb_binary_v2,snp=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/sample_overlap/raw/ukb_binary_v2.qc.snp,script=/mnt/lustre/groups/ukbiobank/usr/sam/project/PRS/sample_overlap/scripts/p2t.R ../scripts/run_gwas.sh
# Need an input file prefix
# Required input
# bfile -> Genotype file prefix
# snp -> SNP filtering file
# file -> PRefix of GWAS guide
# script -> processing script (p2t.Rfor converting binary assoc file
filename=${file}.gwas.${SGE_TASK_ID}
{
    read -r line
    while read -r line; do
        info=( ${line} )
        out=${info[0]}
        pheno_file=${info[3]}
        plink \
            --bfile ${bfile} \
            --extract ${snp} \
            --assoc \
            --autosome \
            --pheno ${pheno_file} \
            --out ${out}
        # Now add the A1 A2 column onto the association file
        # Also compress it to save 5 times the memory space
        awk 'NR==FNR{a[$2]=$5;b[$2]=$6} NR!=FNR && FNR==1{print $0" A1 A2"} NR!=FNR && FNR!=1{print $0,a[$2],b[$2]}' ${bfile}.bim ${out}.qassoc | gzip > ${out}.assoc.gz
        rm  ${out}.qassoc
    done 
} < "$filename"
