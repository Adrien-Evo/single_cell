#!/bin/bash

#$ -l h_vmem=8G
#$ -q default
#$ -cwd
#$ -V
#$ -N clean$JOB_ID
#$ -o clean.out$JOB_ID
#$ -e clean.err$JOB_ID
#$ -m a
#$ -M isabel.alves@oicr.on.ca

#setting the sart of the job
res1=$(date +%s.%N)

folderSC="/.mounts/labs/awadallalab/scratch/SingleCell/vcfs_sc/allSamples"

#loading VCFtools
module load vcftools/0.1.14

tag1=`echo $1 | sed 's/\([^\.]*\)\..*/\1/'`
echo $tag1
#tag2=`echo $2 | sed 's/[^\.]*.\([^\.]*\)\..*/\1/'`
#echo $tag2

vcftools --vcf ${folderSC}/$1 --chr $2 --recode --recode-INFO-all --stdout | \
vcftools --vcf - --exclude-bed ${folderSC}/MappingRegions/blackLists_$2.txt --recode --recode-INFO-all --stdout | \
vcftools --vcf - --bed ${folderSC}/MappingRegions/accessGen_$2.txt --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.$2

if [ $2 == "chr1" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr1p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr1p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr1q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr1q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

	
elif [ $2 == "chr2" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr2p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr2p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr2q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr2q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

	
elif [ $2 == "chr3" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr3p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr3p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr3q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr3q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

elif [ $2 == "chr4" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr4p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr4p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr4q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr4q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q


elif [ $2 == "chr5" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr5p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr5p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr5q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr5q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
elif [ $2 == "chr6" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr6p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr6p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr6q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr6q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q

elif [ $2 == "chr7" ]; 
then
	#getting the start and end of the short arm
	small_start=`grep chr7p ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	small_end=`grep chr7p ${folderSC}/chroms_size_arms | cut -d$'\t' -f3`
	
	#getting the start and end of the longest arm
	large_start=`grep chr7q ${folderSC}/chroms_size_arms | cut -d$'\t' -f2`
	large_end=`grep chr7q ${folderSC}/chroms_size_arms | cut -d$'\t' -f3` 
	
	#dividing vcf file according to the short and longest arm coordinates
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $small_start --to-bp $small_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --chr $2 --from-bp $large_start --to-bp $large_end --recode --recode-INFO-all --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
	#extracting genotypes and QQ
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}p

	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.${2}q
	
else 
	
vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --extract-FORMAT-info GT --out ${folderSC}/$tag1.SingleCellsHetSNPs.$2
vcftools --vcf ${folderSC}/$tag1.SingleCellsHetSNPs.$2.recode.vcf --extract-FORMAT-info GQ --out ${folderSC}/$tag1.SingleCellsHetSNPs.$2

fi

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


