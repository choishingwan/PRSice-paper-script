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
# bfile -> Genotype file prefix
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
loc=`pwd`
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
        if [[ "${herit}" -eq 0 ]]; then 
            # We don't want to analyze data with more than 10000 sample or else it will 
            # blow up our storage space
            if [[ "${target_size}" -lt 100000 ]]; then
                # Need to work within a specific folder to avoid over-writing the results between different jobs
                mkdir -p ${loc}/${out}-tmp
                zcat ${base_assoc} > ${loc}/${out}-tmp/${out}.assoc
                awk 'NR!=1 {print $1,$3}' ${target_pheno} > ${loc}/${out}-tmp/${out}.pheno
                cd ${loc}/${out}-tmp
            
                # Use only 1 thread for PRSice such that it is fair for other programs
                \time -f "%e %S %U %P %K %I %O %W %M" -o ${out}.tmp \
                    Rscript ${prsice} \
                        base ${loc}/${out}-tmp/${out}.assoc \
                        target ${loc}/${target_geno} \
                        covary F \
                        pheno.file ${loc}/${out}-tmp/${out}.pheno \
                        binary.target F \
                        plink /opt/apps/bioinformatics/plink2/1.90b3.38/plink
                # Get R2, 
                r2=`awk -v maxR=0 'NR!=1 && $3>maxR{maxR=$3; maxP=$2}END{print maxR}' PRSice_RAW_RESULTS_DATA.txt`
                pvalue=`awk -v maxR=0 'NR!=1 && $3>maxR{maxR=$3; maxP=$2}END{print maxP}' PRSice_RAW_RESULTS_DATA.txt`
                awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,P,"NA","NA","PRSice-1.25"}' ${out}.tmp > ${loc}/${out}.prsice.old.result
                cd ${loc}
                rm -rf ${loc}/${out}-tmp
            fi
        fi
        
    done } < "$filename"
