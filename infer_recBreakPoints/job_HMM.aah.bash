#!/bin/bash

#$ -l h_vmem=2G
#$ -cwd
#$ -V
#$ -N hmm_$JOB_ID 
#$ -o hmm_$JOB_ID.out
#$ -e hmm_$JOB_ID.err
#$ -m a

#setting the start of the job
res1=$(date +%s.%N)
#module load openmpi
#3 arguments: date, chromosome, and ID

module load R

SCRATCH="/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"
scriptToRun="runningHMM.aah"

chmod +x $SCRATCH/$scriptToRun.R
$SCRATCH/$scriptToRun.R $1 $2 $3

#timing the job
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "dt-86400*$dd" | bc)
dh=$(echo "dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

