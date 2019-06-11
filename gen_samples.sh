#BSUB -L /bin/sh             # script shell language (/bin/tcsh, /bin/ksh etc.)
#BSUB -n 1              # number of cores required
#BSUB -J gen_samples[1-10]               # name of job
#BSUB -q express               # target queue for job execution
#BSUB -W 2:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o %J.stdout
#BSUB -eo %J.stderr


module add R
out=sim${LSB_JOBINDEX}
path=/sc/orga/projects/psychgen/ukb/usrs/sam/project/prs/giga_software/

mkdir -p ${path}/work/${out}
cd ${path}/work/${out}
Rscript ${path}/scripts/paper_master.R -b ${path}/raw/ukb_binary -f ${path}/raw/ukb_binary.qc.fam -s ${path}/raw/ukb_binary.qc.snp -w ${path}/software/BBS/BBS -o ${out} --batch ${path}/raw/ukb1251_cal_chr1_v2_s488363.fam --pc ${path}/raw/ukb1251.pc
