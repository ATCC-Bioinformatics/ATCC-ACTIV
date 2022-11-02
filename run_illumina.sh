#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#paired end data
#for rev in /home/shared/activ/review/illumina/redownload/ERR4439734*_2.fastq; do
#for rev in /home/shared/activ/fastqs/illumina/*_2.fastq; do
for rev in /home/shared/activ/may_2022/illumina/*_2.fastq; do
	prefix=$(basename $rev _2.fastq)
	mkdir -p /home/shared/activ/may_2022/illumina/$prefix/pre_QC
	mkdir -p /home/shared/activ/may_2022/illumina/results/$prefix
#	echo $(dirname $rev)/${prefix}_1.fastq.gz
	bash /home/shared/activ/scripts/april_BEI_illumina.sh $(dirname $rev)/${prefix}_1.fastq $rev /home/shared/activ/may_2022/illumina/results/$prefix/ /home/shared/activ/scripts/NC_045512.2.fasta 48
done
#BEI_illumina_withhost.sh
