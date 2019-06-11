#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 1              # number of cores required
#BSUB -J %FILE%_oprsice[1-80]%36
#BSUB -q express               # target queue for job execution
#BSUB -W 10:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o %FILE%_oprsice.o%J.%I
#BSUB -eo %FILE%_oprsice.e%J.%I
#BSUB -M 40000
module add plink/1.90b6.7
module add R
# Unfortunately, we can't run PRSice-1.25 on the /dev/shm due to the large amount of intermediate file generated
loc=/local/tmp/chois26_prsice${LSB_JOBID}_${LSB_JOBINDEX}

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p ${loc}
chmod 700 ${loc}
cur=`pwd`
prsice=/sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/software/old-prsice/PRSice_v1.25/PRSice_v1.25.R
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
	if [[ $(echo ${herit}'=='0.2 | bc -l) -eq 1 ]]; then
            # We don't want to analyze data with more than 10000 sample or else it will 
            # blow up our storage space
            if [[ "${target_size}" -lt 100000 ]]; then
                # Need to work within a specific folder to avoid over-writing the results between different jobs
                base_name=$(basename -- "${base_assoc}")
        	target_pheno_name=$(basename -- "${target_pheno}")
        	target_geno_name=$(basename -- "${target_geno}")
        	prsice_name=$(basename -- "${prsice}")
        	zcat ${base_assoc} > ${loc}/${base_name}
        	awk 'NR!=1 {print $1,$3}' ${target_pheno} > ${loc}/${target_pheno_name}
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
                cd ${loc}/
        	ls -lht    
                # Use only 1 thread for PRSice such that it is fair for other programs
                \time -f "%e %S %U %P %K %I %O %W %M" -o ${out}.tmp \
                    Rscript ${prsice} \
                        base ${base_assoc} \
                        target ${target_geno} \
                        covary F \
                        pheno.file ${target_pheno} \
                        binary.target F \
                        plink /hpc/packages/minerva-common/plink/1.90b6.7/plink
                # Get R2,
                ls -lht 
                r2=`awk -v maxR=0 'NR!=1 && $3>maxR{maxR=$3; maxP=$2}END{print maxR}' PRSice_RAW_RESULTS_DATA.txt`
                pvalue=`awk -v maxR=0 'NR!=1 && $3>maxR{maxR=$3; maxP=$2}END{print maxP}' PRSice_RAW_RESULTS_DATA.txt`
		awk -v BaseSize=${base_size} -v Size=${target_size} -v Herit=${herit} -v NumCausal=${num_causal} -v R2=${r2} -v P=${pvalue} -v VR2=${valid_r2} -v VP=${valid_pvalue} 'NR==1{print "Real Kernel User Percentage AvgTotalMem Input Output Swap MaxKB Base.Size Target.Size Heritability Num.Causal Target.R2 Target.P Valid.R2 Valid.P Program"}{print $0,BaseSize, Size,Herit,NumCausal,R2,P,VR2,VP,"PRSice-1.25"}' ${out}.tmp > ${out}.prsice.old.result
		cp ${out}.prsice.old.result /sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/work/%FILE%/
		rm -rf ${loc}
            fi
        fi
        
    done } < "$filename"
