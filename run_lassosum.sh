#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 1              # number of cores required
#BSUB -J %FILE%_lassosum[1-80]%30              # name of job
#BSUB -q normal               # target queue for job execution
#BSUB -R "span[hosts=1]"
#BSUB -W 16:00                # wall clock limit for job
#BSUB -M 120000
#BSUB -u choishingwan@gmail.com         # email address for output
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o %FILE%_lassosum.o%J.%I
#BSUB -eo %FILE%_lassosum.e%J.%I
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

module add R
# Need an input file prefix
# Required input
# bfile -> Genotype file prefix
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
lassosum=%LASSOSUM%
LD=%LD%
file=%FILE%
loc=/dev/shm/chois26_lassosum_${LSB_JOBID}_${LSB_JOBINDEX}
mkdir -p ${loc}
chmod 700 ${loc}
# First install lassosum locally to avoid share object error
ls -lh ${loc}
du -h --max-depth=1 ${loc}
filename=${file}.prs.${LSB_JOBINDEX}
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
        lassosum_name=$(basename -- "${lassosum}")
        ld_name=$(basename -- "${LD}")
        cp ${base_assoc} ${loc}/${base_name}
        cp ${target_pheno} ${loc}/${target_pheno_name}
        cp ${lassosum} ${loc}/${lassosum_name}
        cp ${target_geno}.bed ${loc}/${target_geno_name}.bed
        cp ${target_geno}.bim ${loc}/${target_geno_name}.bim
        cp ${target_geno}.fam ${loc}/${target_geno_name}.fam
        cp ${LD} ${loc}/${ld_name}
        base_assoc=${loc}/${base_name}
        target_pheno=${loc}/${target_pheno_name}
        target_geno=${loc}/${target_geno_name}
        final_out=${out}
        LD=${loc}/${ld_name}
        lassosum=${loc}/${lassosum_name}
        out=${loc}/${out}
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
		echo ${r2}
		ls ${out}*
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
                awk -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v VR2=${valid_r2}  'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize,Size,Herit,NumCausal,R2,"NA",VR2,"NA","lassosum"}' ${out}.lassosum.tmp > ${final_out}.lassosum.result
            else
				echo ${r2} 
                awk  -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2}   'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,"NA","NA","NA","lassosum"}' ${out}.lassosum.tmp > ${final_out}.lassosum.result
            fi
            rm ${out}.Rdata
            
            rm ${out}.lassosum.tmp
        fi
    done } < "$filename"
