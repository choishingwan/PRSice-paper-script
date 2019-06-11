#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 1              # number of cores required
#BSUB -J gen_plink[1-20]               # name of job
#BSUB -q express               # target queue for job execution
#BSUB -W 5:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

module add plink/1.90b6.7
# Need an input file prefix
# Required input
# bfile -> Genotype file prefix
# snp -> SNP filtering file
# file -> PRefix of GWAS guide
# script -> processing script (p2t.Rfor converting binary assoc file
filename=%FILE%.plink.${LSB_JOBINDEX}
{
    read -r line
    while read -r line; do
        info=( ${line} )
        out=${info[0]}
        plink \
            --bfile /sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/raw/ukb_binary \
            --extract /sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/raw/ukb_binary.qc.snp \
            --pheno /sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/raw/ukb_binary.ldpred.fam \
	    --keep ${out} \
            --make-bed \
            --out ${out}
        # Now add the A1 A2 column onto the association file
        # Also compress it to save 5 times the memory space
    done 
} < "$filename"
