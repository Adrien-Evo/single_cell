#!/bin/bash

chrID=$1
date="Dec27"

while IFS='' read -r line || [[ -n "$line" ]];
do

indID=`echo $line | cut -d$' ' -f2 | cut -d$'.' -f2`;
echo $indID;

#launching phasing
qsub -pe smp 10 /.mounts/labs/awadallalab/scratch/ialves/wgs_sc/job_HQ25_10_5.bash $indID $date $chrID

done < "$HOME/files_Dec17/fileNames_SC.txt"
