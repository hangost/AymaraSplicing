# Mapping for RNA-seq data

#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o ./tmp/STAR.o.%j
#SBATCH --err ./tmp/STAR.e.%j

ml R/3.3.2
ml samtools
ml star/2.5
sampleNM=${1##*/}
sampleNM=${sampleNM%%_*}
Foutdir=${4}/${sampleNM}/

echo ${Foutdir}
mkdir -p ${Foutdir}

STAR --genomeDir ./STAR-genome2/STAR/Total/ --sjdbFileChrStartEnd ./STAR-genome2/STAR/Total/sjdbList.fromGTF.out.tab --runThreadN $3 --outFileNamePrefix ${Foutdir} --outSAMtype BAM SortedByCoordinate --readFilesIn $1 $2
