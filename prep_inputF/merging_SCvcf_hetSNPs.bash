#!/bin/bash

#$ -l h_vmem=8G
#$ -q default
#$ -cwd
#$ -V
#$ -N vcf$JOB_ID
#$ -o vcf.out$JOB_ID
#$ -e vcf.err$JOB_ID
#$ -m a
#$ -M isabel.alves@oicr.on.ca

#setting the sart of the job
res1=$(date +%s.%N)

folderSC="/.mounts/labs/awadallalab/scratch/SingleCell/vcfs_sc/allSamples"
folderWGS="/.mounts/labs/awadallalab/scratch/SingleCell/WGS_BULK_common.hetSNPs"

#loading VCFtools
module load vcftools/0.1.14
module load bcftools-1.7/1.7

tag1=`echo $1 | sed 's/\([^\.]*\)\..*/\1/'`
echo $tag1
tag2=`echo $2 | sed 's/[^\.]*.\([^\.]*\)\..*/\1/'`
echo $tag2

# if [ ! -d "${folderWGS}/highM_vcf" ]; 
# 	then
# 		mkdir ${folderWGS}/highM_vcf
# fi
# 
# headerSize=`grep -n "#CHROM" ${folderWGS}/$2 | cut -d$':' -f1` 
# 
# sed -n 1,${headerSize}p ${folderWGS}/$2 > ${folderWGS}/highM_vcf/$tag2.HM.vcf
# 
# grep "MAPPABILITY=1" ${folderWGS}/$2 >> ${folderWGS}/highM_vcf/$tag2.HM.vcf

vcftools --gzvcf ${folderSC}/$1 --gzdiff ${folderWGS}/$2 --diff-site --out ${folderSC}/$tag2.hetSNPs.SC

grep "B" ${folderSC}/$tag2.hetSNPs.SC.diff.sites_in_files | cut -d$'\t' -f1-2 > ${folderSC}/$tag2.hetSNPs.SC.InSingleCells.txt

vcftools --gzvcf ${folderSC}/$1 --positions ${folderSC}/$tag2.hetSNPs.SC.InSingleCells.txt --recode --recode-INFO-all --out ${folderSC}/$tag2.SingleCellsHetSNPs

#rm ${folderSC}/$tag2.hetSNPs.HM.SC.diff.sites_in_files
#rm ${folderSC}/$tag2.hetSNPs.HM.InSingleCells.txt

#timing the job
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

echo "Total runtime: $dd days $dh hrs $dm min $ds secs"


