#!/bin/bash

#$ -l h_vmem=4G
#$ -cwd
#$ -V
#$ -N jobID.$JOB_ID
#$ -o jobID.$JOB_ID.out
#$ -e jobID.$JOB_ID.err
#$ -m a
#$ -M isabel.alves@oicr.on.ca

#setting the start of the job
res1=$(date +%s.%N)
#module load openmpi
module load R

SCRATCH="/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"
scriptToRun="visualizeBreakPoints"

# date; chr; idID

chmod +x $SCRATCH/$scriptToRun.R
Rscript --vanilla $SCRATCH/$scriptToRun.R $1 $2 $3

#timing the job
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "dt-86400*$dd" | bc)
dh=$(echo "dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

