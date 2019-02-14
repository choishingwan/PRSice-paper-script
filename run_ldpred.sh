#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q HighMemLongterm.q,LowMemLongterm.q,LowMemShortterm.q
#$ -t 1:80
#$ -l h_rt=72:00:00
#$ -l h_vmem=100G


module add bioinformatics/plink2/1.90b3.38
module add general/python/3.5.1
source ~/.bash_profile 
# Need an input file prefix
# Required input
# bfile -> Genotype file prefix
# file -> Prefix of PRS guide
# ldpred -> location of LDPRed

# 0 = Base.Assoc
# 1 = Base.Size 
# 2 = Prefix
# 3 = Target.Pheno
# 4 = Target.Size
# 5 = Validate.Pheno 
# 6 = Num.Causal 
# 7 = Heritability 
# 8 = Radius 
# 9 = Target genotype file
# 10 = Valid genotype file
filename=${file}.prs.${SGE_TASK_ID}
{
    read -r line
    while read -r line; do
        info=( ${line} )
        base_assoc=${info[0]}
        base_size=${info[1]}
        out=${info[2]}.ldpred
        target_pheno=${info[3]}
        valid_pheno=${info[5]}
        target_size=${info[4]}
        num_causal=${info[6]}
        herit=${info[7]}
        radius=${info[8]}
        target_geno=${info[9]}
        valid_geno=${info[10]}
		rm ${out}.coord
        # Ignore heritability = 0 cases as those will crash LDpred + one shouldn't
        # perform PRS analysis on traits with SNP heritability = 0
        if [[ "${herit}" -ne 0 ]]; then
            \time -f "%e %S %U %P %K %I %O %W %M" -o ${out}.ldpred.tmp \
            sh -c "
                python3.5 ${ldpred} coord \
                    --rs SNP \
                    --A1 A1 \
                    --A2 A2 \
                    --pos BP \
                    --chr CHR \
                    --pval P \
                    --eff BETA \
                    --beta \
                    --ssf-format CUSTOM \
                    --N ${base_size} \
                    --ssf ${base_assoc} \
                    --out ${out}.coord \
                    --gf ${target_geno};

                python3.5 ${ldpred} gibbs \
                    --cf ${out}.coord \
                    --ldr ${radius} \
                    --ldf ${out}.ld \
                    --out ${out}.weight \
                    --N ${base_size};

                python3.5 ${ldpred} score \
                    --gf ${target_geno} \
                    --rf ${out}.weight \
                    --out ${out}.score \
                    --pf ${target_pheno} \
                    --pf-format LSTANDARD 
    "
            result=(`awk -v max=-1,id=0 '$2>max{max=$2;id=$1}END{print id,max}' ${out}.score`)
            if [[ "${result[0]}" == "inf" ]]; then
                python3.5 ${ldpred} inf \
                    --cf ${out}.coord \
                    --ldr ${radius} \
                    --ldf ${out}.valid.ld \
                    --out ${out}.weight \
                    --N ${base_size}
                    ##_LDpred-inf.txt
                rm ${out}.coord
                python3.5 ${ldpred} score \
                    --gf ${valid_geno} \
                    --rf ${out}.weight \
                    --out ${out}.score \
                    --pf ${valid_pheno} \
                    --pf-format LSTANDARD  
            else
                python3.5 ${ldpred} gibbs \
                    --cf ${out}.coord \
                    --ldr ${radius} \
                    --ldf ${out}.valid.ld \
                    --out ${out}.weight \
                    --N ${base_size} \
                    --f ${result[0]}
                rm ${out}.ldpred_LDpred-inf.txt 
                rm ${out}.coord
                python3.5 ${ldpred} score \
                    --gf ${valid_geno} \
                    --rf ${out}.weight \
                    --out ${out}.score \
                    --pf ${valid_pheno} \
                    --pf-format LSTANDARD 
            fi
            vres=(`awk -v max=-1,id=0 '$2>max{max=$2;id=$1}END{print id,max}' ${out}.vldpred`)
            awk -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${result[1]} -v VR2=${vres[1]} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,Size,Herit,NumCausal,R2,"NA",VR2,"NA","LDPred"}' ${out}.ldpred.tmp > ${out}.ldpred.result
            rm ${out}.ldpred.tmp
        fi
    done } < "$filename"
