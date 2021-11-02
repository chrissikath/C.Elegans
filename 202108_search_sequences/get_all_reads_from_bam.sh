#!/bin/sh
# count all reads in given bam file
FILE="/home/udo/worms/*.bam"
set -e -x

touch "readsinBam.csv"
for i in $FILE; do
	echo $i
	total_counts=$(samtools view -c $i) #count reads 
	echo "${i##*/}\t $total_counts" >> readsinBam.csv
done
