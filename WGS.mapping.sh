#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o ./tmp/Hisat.o.%j
#SBATCH --err ./tmp/Hisat.e.%j

ml R/3.3.2
ml samtools
ml star/2.5
ml hisat2/2.0.4


for_read=$1
rev_read=${1%%1.txt}2.txt

sampleNM=${1##*/}
sampleNM=${sampleNM%%_*}
Foutdir=./9.WGS.Hisat2/${sampleNM}/

echo $sampleNM
mkdir -p ${Foutdir}

hisat2 -p 12 -x ../HISAT2-genome/hg19.hisat2 -1 $for_read -2 $rev_read -S ${Foutdir}/${sampleNM}.sam
