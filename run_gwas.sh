#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 1              # number of cores required
#BSUB -J run_gwas[1-20]               # name of job
#BSUB -q express               # target queue for job execution
#BSUB -W 2:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

module add plink/1.90b6.7
file=%FILE%
# Need an input file prefix
# Required input
# bfile -> Genotype file prefix
# snp -> SNP filtering file
# file -> PRefix of GWAS guide
# script -> processing script (p2t.Rfor converting binary assoc file
bfile=/sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/raw/ukb_binary 
snp=/sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/raw/ukb_binary.qc.snp

filename=${file}.gwas.${LSB_JOBINDEX}
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
