#!/usr/bin/bash
#SBATCH --nodelist=n003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#paired end data
#for rev in /home/shared/activ/fastqs/ont/*.fastq.gz; do
for rev in /home/shared/activ/fastqs/ont/*gz; do
	prefix=$(basename $rev .fastq.gz)
	mkdir -p /home/shared/activ/zresults_nanopore_medaka/$prefix/pre_QC
#	echo $(dirname $rev)/${prefix}_1.fastq.gz
	bash zBEI_nanopore_medakaNanocall.sh $rev /home/shared/activ/zresults_nanopore_medaka/$prefix /home/shared/activ/scripts/NC_045512.2.fasta 48
done
