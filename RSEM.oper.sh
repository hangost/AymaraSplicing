#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o ./tmp/rsem.o.%j
#SBATCH --err ./tmp/rsem.e.%j

ml R/3.3.2
ml samtools
ml star/2.5
ml rsem/1.3.0

sampleNM=${1##*/}
sampleNM=${sampleNM%%_*}
Foutdir=${4}/${sampleNM}

echo ${Foutdir}
mkdir -p ${Foutdir}

rsem-calculate-expression --paired-end -p $3 --star $1 $2 ./4.rsem/WithOutpatch/rsem.hg19 ${Foutdir}/${sampleNM}.result
