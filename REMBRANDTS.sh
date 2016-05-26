#!/bin/bash

####################### define executables
deseq="./src/_R/DESeq.R"
rambrandts="./src/_R/REMBRANDTS.R"

####################### identify the input arguments
jobid=$1
metadata=$2
refdir=$3
stringency=$4

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash MoSBAT.sh <jobID> <motif1.pwm> <motif2.pwm> <motif_type> <sequence_length> <sequence_count>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input metadata file: "$metadata
echo "Reference directory for HTSeq-Count files: "$refdir
echo "Stringency for filtering measurements: "$stringency

if [ -e "$metadata" ]; then
	echo "Metadata file found."
else
	echo "ERROR: Metadata file was not found."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder

####################### run DESeq
echo -n -e "GeneID\t" > $tmp_folder/"vsd_normalized.exonic.all.centered.mx.txt"
echo -n -e "GeneID\t" > $tmp_folder/"vsd_normalized.exonic.all.mx.txt"
echo -n -e "GeneID\t" > $tmp_folder/"vsd_normalized.intronic.all.centered.mx.txt"
echo -n -e "GeneID\t" > $tmp_folder/"vsd_normalized.intronic.all.mx.txt"

Rscript $deseq $jobid $metadata $refdir


####################### define output path
out_folder="./out/"$jobid
mkdir -p $out_folder

####################### run REMBRANDTS

Rscript $rambrandts $jobid $stringency