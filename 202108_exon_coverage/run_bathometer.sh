#!/bin/sh
set -e -x

#test worked:
#bathometer index G1-N2-.bam.idx /home/udo/worms/G1-N2-.bam 

#make indexes for all bam files
#for i in /home/udo/worms/*.bam; do
#	echo ${i##*/}.idx
#	bathometer index ${i##*/}.idx $i
#ydone 
#worked!

#count exon coverage for all exons in bed file
 
for i in *.bed; do
	echo $i
	bathometer fpkm *.idx ::: $i > exon_coverage_${i%.bed}_FPKM.tsv
done

