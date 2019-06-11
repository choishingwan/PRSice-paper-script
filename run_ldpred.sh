#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 2              # number of cores required
#BSUB -J %FILE%_ldpred[1-80]%20               # name of job
#BSUB -R "span[hosts=1]"
#BSUB -q normal               # target queue for job execution
#BSUB -W 24:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o run_ldpred.o%J.%I
#BSUB -eo run_ldpred.e%J.%I
#BSUB -M 110000000
module add python/3.5.0
ldpred=%LDPRED%
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
loc=/dev/shm/chois26_ldpred${LSB_JOBID}_${LSB_JOBINDEX}
#loc=/local/tmp/chois26_ldpred${LSB_JOBID}_${LSB_JOBINDEX}
mkdir -p ${loc}
chmod 700 ${loc}
ls -lht ${loc}
du -h --max-depth=1 ${loc}
filename=%FILE%.prs.${LSB_JOBINDEX}
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
        radius=${info[8]}
        target_geno=${info[9]}
        valid_geno=${info[10]}
        base_name=$(basename -- "${base_assoc}")
        target_pheno_name=$(basename -- "${target_pheno}")
        target_geno_name=$(basename -- "${target_geno}")
        ldpred_name=$(basename -- "${ldpred}")
        ldpred_dir=$(dirname "${ldpred}")
        cp ${base_assoc} ${loc}/${base_name}
        cp ${target_pheno} ${loc}/${target_pheno_name}
        cp -r ${ldpred_dir} ${loc}
        cp ${target_geno}.bed ${loc}/${target_geno_name}.bed
        cp ${target_geno}.bim ${loc}/${target_geno_name}.bim
        cp ${target_geno}.fam ${loc}/${target_geno_name}.fam
        base_assoc=${loc}/${base_name}
        target_pheno=${loc}/${target_pheno_name}
        target_geno=${loc}/${target_geno_name}
        final_out=${out}
        ldpred=${loc}/ldpred/LDpred.py
        out=${loc}/${out}
		rm ${out}.coord
        # Ignore heritability = 0 cases as those will crash LDpred + one shouldn't
        # perform PRS analysis on traits with SNP heritability = 0
        if (( $(echo ${herit}'=='0 | bc -l) )); then 
            echo "Skip Heritability == 0"
        else
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
                    --pf-format LSTANDARD >${out}.log 2>&1 
    "
			result=`grep highest ${out}.log | grep -v Method | awk '{print $7}' | sed 's/,//g'`
            #result=(`awk -v max=-1,id=0 '$2>max{max=$2;id=$1}END{print id,max}' ${out}.score`)
            if [[ "${valid_geno}" != "NA" ]]; then
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
                awk -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${result} -v VR2=${vres[1]} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,"NA",VR2,"NA","LDpred"}' ${out}.ldpred.tmp > ${final_out}.ldpred.result
            else 
                awk -v Size=${target_size} -v BaseSize=${base_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${result} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,"NA","NA","NA","LDpred"}' ${out}.ldpred.tmp > ${final_out}.ldpred.result
            fi
            rm ${out}.ldpred.tmp ${out}.*weight* ${out}.*score*
        fi
    done } < "$filename"

rm -rf ${loc}
