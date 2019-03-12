#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q HighMemLongterm.q,LowMemLongterm.q,LowMemShortterm.q
#$ -t 1:80
#$ -l h_rt=72:00:00

module add bioinformatics/plink2/1.90b3.38
module add bioinformatics/R/3.3.3
module add compilers/gcc/6.2.0
# Need an input file prefix
# Required input
# file -> Prefix of PRS guide
# PRSice -> location of PRSice

# 0 = Base.Assoc
# 1 = Base.Size 
# 2 = Prefix
# 3 = Target.Pheno
# 4 = Target.Size
# 5 = Validate.Pheno 
# 6 = Num.Causal 
# 7 = Heritability 
# 8 = Radius 
# 9 = Target Genotype
# 10 = Valid Genotype

filename=${file}.prs.${SGE_TASK_ID}
{
    read -r line
    while read -r line; do
        info=( ${line} )
        base_assoc=${info[0]}
        out=${info[2]}
        target_pheno=${info[3]}
        valid_pheno=${info[5]}
        target_size=${info[4]}
        num_causal=${info[6]}
        herit=${info[7]}
        target_geno=${info[9]}
        valid_geno=${info[10]}
        # Use only 1 thread for PRSice such that it is fair for other programs
        \time -f "%e %S %U %P %K %I %O %W %M" -o ${out}.tmp \
            ${prsice} \
                --base ${base_assoc} \
                --target ${target_geno} \
                --pheno-file ${target_pheno} \
                --thread 1 \
                --binary-target F \
                --out ${out} \
                --stat T \
                --beta
        rm ${out}.best
        rm ${out}.prsice
        threshold=`awk 'NR!=1{print $3}' ${out}.summary`
        r2=`awk 'NR!=1{print $4}' ${out}.summary`
        pvalue=`awk 'NR!=1{print $10}' ${out}.summary`
        if [[ "${valid_geno}" != "NA" ]]; then
            ${prsice} \
                --base ${base_assoc} \
                --target ${valid_geno} \
                --binary-target F \
                --bar-levels ${threshold} \
                --no-full \
                --beta \
                --fastscore \
                --out ${out}.valid \
                --pheno-file ${valid_pheno} \
                --thread max
            rm ${out}.valid.best
            rm ${out}.valid.prsice
            valid_r2=`awk 'NR!=1{print $4}' ${out}.valid.summary`
            valid_pvalue=`awk 'NR!=1{print $10}' ${out}.valid.summary`
            awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} -v VR2=${valid_r2} -v VP=${valid_pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,P,VR2,VP,"PRSice-2"}' ${out}.tmp > ${out}.prsice.result
        else
            awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,P,"NA","NA","PRSice-2"}' ${out}.tmp > ${out}.prsice.result
        fi
        rm ${out}.tmp
    done } < "$filename"
