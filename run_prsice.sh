#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 6              # number of cores required
#BSUB -J run_prsice[1-80]%36               # Limit to 36, so we only use 3 host at the same time
#BSUB -R "span[hosts=1]"
#BSUB -q express               # target queue for job execution
#BSUB -W 2:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o run_prsice.o%J.%I
#BSUB -eo run_prsice.e%J.%I
#BSUB -M 30000
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
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
prsice=%PRSICE%
loc=/dev/shm/chois26_prsice_${LSB_JOBID}_${LSB_JOBINDEX}
mkdir -p ${loc}
chmod 700 ${loc}
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
        target_geno=${info[9]}
        valid_geno=${info[10]}
        base_name=$(basename -- "${base_assoc}")
        target_pheno_name=$(basename -- "${target_pheno}")
        target_geno_name=$(basename -- "${target_geno}")
        prsice_name=$(basename -- "${prsice}")
        cp ${base_assoc} ${loc}/${base_name}
        cp ${target_pheno} ${loc}/${target_pheno_name}
        cp ${prsice} ${loc}/${prsice_name}
        cp ${target_geno}.bed ${loc}/${target_geno_name}.bed
        cp ${target_geno}.bim ${loc}/${target_geno_name}.bim
        cp ${target_geno}.fam ${loc}/${target_geno_name}.fam
        prsice=${loc}/${prsice_name}
        base_assoc=${loc}/${base_name}
        target_pheno=${loc}/${target_pheno_name}
        target_geno=${loc}/${target_geno_name}
        final_out=${out}
        out=${loc}/${out}
        ls -lht ${loc}
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
            awk -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} -v VR2=${valid_r2} -v VP=${valid_pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,P,VR2,VP,"PRSice-2"}' ${out}.tmp > ${final_out}.prsice.result
        else
            awk -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,P,"NA","NA","PRSice-2"}' ${out}.tmp > ${final_out}.prsice.result
        fi
        rm ${out}.tmp
        rm ${out}.summary
	cp ${out}.log ${final_out}.log
    done } < "$filename"

rm -rf ${loc}
