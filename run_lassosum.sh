#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q HighMemLongterm.q,LowMemLongterm.q,LowMemShortterm.q
#$ -t 1:80
#$ -l h_rt=72:00:00
#$ -l h_vmem=20G

module add bioinformatics/plink2/1.90b3.38
module add bioinformatics/R/3.3.3
module add compilers/gcc/6.2.0
# Need an input file prefix
# Required input
# file -> Prefix of PRS guide
# lassosum -> location of lassosum.R
# LD -> location of the LD file

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
        base_size=${info[1]}
        out=${info[2]}
        target_pheno=${info[3]}
        valid_pheno=${info[5]}
        target_size=${info[4]}
        num_causal=${info[6]}
        herit=${info[7]}
        target_geno=${info[9]}
        valid_geno=${info[10]}

        if (( $(echo ${herit}'=='0 | bc -l) )); then 
            echo "Skip Heritability == 0"
        else
            # We do not time the validation step
            \time -f "%e %S %U %P %K %I %O %W %M" -o ${out}.lassosum.tmp \
                Rscript ${lassosum} \
                    --base ${base_assoc} \
                    --ref ${target_geno} \
                    --target ${target_pheno} \
                    --genotype ${target_geno} \
                    --ld ${LD} \
                    --size ${base_size} \
                    --out ${out}
                
            r2=`awk 'NR!=1{print $2}' ${out}.lassosum`
            if [[ "${valid_geno}" != "NA" ]]; then
                Rscript ${lassosum} \
                    --base ${base_assoc} \
                    --ref ${valid_geno} \
                    --target ${out}.Rdata \
                    --valid ${valid_pheno} \
                    --genotype ${valid_geno} \
                    --ld ${LD} \
                    --size ${base_size} \
                    --out ${out}
                valid_r2=`awk 'NR!=1{print $3}' ${out}.lassosum`
                awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v VR2=${valid_r2}  'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,"NA",VR2,"NA","lassosum"}' ${out}.lassosum.tmp > ${out}.lassosum.result
            else
				echo ${r2} 
                awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2}   'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,"NA","NA","NA","lassosum"}' ${out}.lassosum.tmp > ${out}.lassosum.result
            fi
            rm ${out}.Rdata
            
            rm ${out}.lassosum.tmp
        fi
    done } < "$filename"
