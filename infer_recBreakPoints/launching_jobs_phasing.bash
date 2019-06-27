#!/bin/bash

#the name of the bash script to be copied 
#and changed has to be passed as an argument.
scriptName=$1
indName="C133"
QCTag="HQ25_16_8"
SCRATCH="/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"
rScriptName="phasing_${QCTag}_May2"
todayDate="May2"

for i in `seq 11 11`; #copies and changes the $scriptName file as many times as there are chromosomes 

	do
		chrNb="chr$i"
		echo "Launching job for ind: $indName chr: $chrNb"
		
		
		cp $SCRATCH/$scriptName.bash $SCRATCH/$scriptName.$indName.$chrNb.bash
		cp $SCRATCH/$rScriptName.R $SCRATCH/$rScriptName.$indName.$chrNb.R
		
		#replace pattern CHRID by the variable $chrNb within the new copied bash file
		sed -i s/JOBID/${indName}_${chrNb}/ $SCRATCH/$scriptName.$indName.$chrNb.bash
		sed -i s/SCRIPTTORUN/${rScriptName}.${indName}.${chrNb}/ $SCRATCH/$scriptName.$indName.$chrNb.bash

		#replace pattern CHRID by the variable $chrNb within the new copied bash file
		sed -i s/TODAY/${todayDate}/ $SCRATCH/$rScriptName.$indName.$chrNb.R
		sed -i s/INDIVIDUAL/${indName}/ $SCRATCH/$rScriptName.$indName.$chrNb.R
		sed -i s/CHROM/${chrNb}/ $SCRATCH/$rScriptName.$indName.$chrNb.R

						
		#make the bash file exectutable
		chmod +x $SCRATCH/$scriptName.$indName.$chrNb.bash
		#runs new chromosome-specific bash file
		#qsub -pe smp 10 $SCRATCH/$scriptName.$indName.$chrNb.bash
		
	done
	
