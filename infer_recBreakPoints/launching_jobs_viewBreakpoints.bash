#!/bin/bash

#the name of the bash script to be copied 
#and changed has to be passed as an argument.
scriptName=$1
rScriptNames="visualizeBreakPoints"
SCRATCH="/.mounts/labs/awadallalab/scratch/ialves/wgs_sc"
dateToday="Mar15"
indID="C24"

for i in `seq 11 11`; #copies and changes the $scriptName file as many times as there are chromosomes 

	do
		chrNb="chr$i"
		echo "Launching job for chr: $chrNb"
		
		cp $SCRATCH/$rScriptNames.R $SCRATCH/$rScriptNames.$indID.$chrNb.R
		cp $SCRATCH/$scriptName.bash $SCRATCH/$scriptName.$indID.$chrNb.bash
		
		#replace pattern CHRID by the variable $chrNb within the new copied bash file
		sed -i s/TODAY/$dateToday/ $SCRATCH/$rScriptNames.$indID.$chrNb.R
		sed -i s/CHROM/$chrNb/ $SCRATCH/$rScriptNames.$indID.$chrNb.R
		sed -i s/INDIVIDUAL/$indID/ $SCRATCH/$rScriptNames.$indID.$chrNb.R
		
		sed -i s/jobID/graphs_${chrNb}/ $SCRATCH/$scriptName.$indID.$chrNb.bash
		sed -i s/SCRIPTTORUN/$rScriptNames.$indID.$chrNb/ $SCRATCH/$scriptName.$indID.$chrNb.bash
		#make the bash file exectutable
		chmod +x $SCRATCH/$scriptName.$indID.$chrNb.bash
		#runs new chromosome-specific bash file
		#qsub $SCRATCH/$scriptName.$indID.$chrNb.bash
		
	done
	
