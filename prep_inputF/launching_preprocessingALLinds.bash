#!/bin/bash

chrID=$1
date="Dec27"


while IFS='' read -r line || [[ -n "$line" ]];
do

indID=`echo $line | cut -d$' ' -f2 | cut -d$'.' -f2`;
echo $indID;
#launching cleaning of the SC files
#qsub merging_SCvcf_hetSNPs.bash $line

#launching subsetting per chromosome
sc_file=`echo $line | cut -d$' ' -f1`
qsub -hold_jid "vcf*" subsetting_chr_extracting.bash $indID.SingleCellsHetSNPs.recode.vcf $chrID

done < "$HOME/files_Dec17/fileNames_SC.txt"
