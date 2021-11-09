#!/bin/sh

#Author Christina Kuhn @github: https://github.com/chrissikath
# pre-used scripts from Alex 
# Shell Script, welches anhand von RNA-Seq Dateien Splice-Varianten findet Part 2
# R version: 3.6
# BiocManager::install(version = "3.10")
# R-Packages installieren 

#1.0 files
GENOME="/home/horn/2021_Celegans/20210902_alex_data/genome/wormbase/caenorhabditis_elegans.PRJNA13758.WBPS10.genomic.fa"
exon_information_script="/home/christina/C_elegans/analysis/20211022_EasySpliceVariantAnalyses/exon-information.r"
packages_functions_script="/home/christina/C_elegans/analysis/20211022_EasySpliceVariantAnalyses/packages_functions_Alex.R"

#2.0 files
motifs_gene_file="motifs_adgrl1"
set_working_dir="/home/christina/C_elegans/analysis/20211022_EasySpliceVariantAnalyses/"
visualization_splice_variants_script="/home/christina/C_elegans/analysis/20211022_EasySpliceVariantAnalyses/visualization_splice_variants.R" 

#Gene information
gene="Lat-1" # define gene name
gene_name="ADGRL1" 
start="8896841" # define start
stop="8908666"   # define stop
chr="II" #define chromosome
strand="+" #define strand 

############################## 1.0 shell part ###################

# exon-based comparison
cd  results
R --vanilla < $exon_information_script $gene_name $packages_functions_script \
      `ls | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}'`
cd ..

# concat all appended files
cat results/*_append.gtf > all_transcripts.gtf
# make bed file from exons.gtfs
cat results/unique_exons.gtf | awk 'OFS="\t" {print $1,$4-1,$5,$12,$7}' | tr -d '";' > results/unique_exons.bed

bedtools getfasta -fi ${GENOME} \
	-bed results/unique_exons.bed \
	-fo exon_seqs.fa -name

# get sequence of whole locus 
sed -e 1b -e '$!d' results/unique_exons.bed > locus.tmp
first=`cat locus.tmp | head -1 | awk -v x=2 '{print $x}'`
last=`cat locus.tmp | sed -n 2p | awk -v x=3 '{print $x}'`
rm locus.tmp

echo $chromosome $first $last > whole_locus.tmp
cat whole_locus.tmp | awk -v OFS='\t' '{print $1, $2,$3}' > whole_locus.bed

rm whole_locus.tmp

bedtools getfasta -fi ${GENOME} \
	-bed whole_locus.bed \
	-fo whole_locus.fa -name
	
############################# 2.0 R Part ##########################################

# Wichtige Output Dateien:
	# transcripts variants pdf 
	# sequence_and_orf_analysis pdf 
	# orfs.fa
R  --vanilla < $visualization_splice_variants_script $gene_name $motifs_gene_file $packages_functions_script $set_working_dir