#!/bin/sh
set -e -x	# Aussteigen beim 1.ten Fehler und Fehler ausgeben
#$ chmod a+rx my-script.sh
#$ ./my-script.sh

# get depth at regions .bed, these are FPKM
# bathometer fpkm [index.idx...] ::: [region|regions.bed|regions.gff]

# first do the indexing
#---------------------------------

#for i in SRR7961322*.bam; do
#    bathometer index ${i%.bam}.idx $i
#done

# works


# do it for all .bam
#---------------------------------

#for i in *.bam; do
#    bathometer index ${i%.bam}.idx $i
#done


# then do query
#---------------------------------

for i in *.bed; do
    bathometer fpkm /home/horn/2019_Suntsova_normal_tissues/hisat2/*.idx :::  $i >Suntsova_2019_${i%.bed}_FPKM.tsv
done
