#!/bin/sh

#Author Christina Kuhn @github: https://github.com/chrissikath
# pre-used scripts from Alex 
# Shell Script, welches anhand von RNA-Seq Dateien Splice-Varianten findet Part 1
# STAR version 2.7.4a
# Stringtie version 1.3.3b

GENOME="/home/horn/2021_Celegans/20210902_alex_data/genome/wormbase/caenorhabditis_elegans.PRJNA13758.WBPS10.genomic.fa" 
ANNO_GTF="/home/horn/2021_Celegans/20210902_alex_data/genome/wormbase/caenorhabditis_elegans.PRJNA13758.WBPS10.canonical_geneset.gtf"
STARINDEX_OUTPUT="/home/christina/C_elegans/analysis/20211022_EasySpliceVariantAnalyses/STARIndex" #wo der StarIndex hinsoll
RAW_SEQ="/home/biochemistry/2003KNO-0044/01.RawData/" #Ordner wo raw sequences liegen
SAMPLE_LIST="sample_list_nur_ein_sample.txt" #File welche Dateien benutzt werden sollen 

gene="Lat-1" # define gene name
start="8896841" # define start
stop="8908666"   # define stop
chr="II" #define chromosome
strand="+" #define strand 

############################# 1.1 STAR INDEX ################################

# creates STAR index into directory ${STARINDEX_OUTPUT}
echo '...............create STAR index'; 
STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir ${STARINDEX_OUTPUT} \
-- genomeFastaFiles ${GENOME} \
-- sjdbGTFfile ${ANNO_GTF} \
-- sjdbOverhang 100 -- genomeSAindexNbases 12;

 
##############################################################################

############################ 1.2 Mapping #####################################

# Mapping a list FASTA-files to the STARindex ${STARINDEX_OUTPUT} 
# raw data are zipped fasta files and paired-end reads
# The sample_list file must be created beforehand (first column = name of the FASTA-file 
# without the "_1/2", second column = arbitrary identifier, for example "d1_adult_whole_1") ...
# Example: 	G1-N2	d1_adult_gonad_wildtype_1
#			G2-N2	d1_adult_gonad_wildtype_2
#			G3-N2	d1_adult_gonad_wildtype_3
#			G4-N2	d1_adult_gonad_wildtype_4
#			G5-N2	d1_adult_gonad_wildtype_5

echo '...............mapping: ' ${SAMPLE_LIST};
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample `date +%T`;
	sra=`grep $sample ${SAMPLE_LIST} | cut -f1`;

    #location of the raw data files must be defined
	forw=`ls ${RAW_SEQ}${sra}_1.fastq.gz`;
	rev=`ls ${RAW_SEQ}${sra}_2.fastq.gz`;
	echo $forw;
	echo $rev;
		  
	#sub-directories for the resultinig files of each sample are created (only needed in first run)
	mkdir -p star/$sample/;
	STAR itself
	STAR \
	--runThreadN 20 \
	--genomeDir ${STARINDEX_OUTPUT} \
	--readFilesCommand zcat \
	--readFilesIn $forw $rev \
	--outFileNamePrefix ${PWD}/star/$sample/ \
	--outSAMtype BAM Unsorted \
	--quantMode GeneCounts \
	--outSAMstrandField intronMotif
done;
    
################################################################################

############################ 1.3 Sort bed files ######################
echo '...............sort';
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample;
	samtools sort ${PWD}/star/$sample/Aligned.out.bam \
		-o ${PWD}/star/$sample/Aligned.out.bam_sorted.bam;
	done; 
################################################################################

############################ 2.0 Stingtie ######################################
echo '...............stringtie';
mkdir -p stringtie/
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample `date +%T`;
	mkdir ${PWD}/stringtie/$sample/;
	stringtie -c 0.1 -f 0.0 -m 50 -a 1 -p 20 ${PWD}/star/$sample/Aligned.out.bam_sorted.bam \
	  -o ${PWD}/stringtie/$sample/$sample.gtf;
	# -c -> minimum read coverage 
	# -f -> disable trimming at the end of the transcripts 
	# -m -> minimum length allowed for the predicted transcripts (normal=200)
	# -a-> filtered out junctions that dont have spliced reads align with at least 1
	# -p -> number of processing threads 
done; 

#################################################################################

############################ 3.0 Extract Gene from gtf results ##################
echo  '...............extract gene from gtf';
mkdir -p results/
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample 'extract gene from gtf' `date +%T` ;
	cat ${PWD}/stringtie/$sample/$sample.gtf \
	 | tail -n+3 | grep ${strand} | awk -v var="${chr}" '{if ($1 == var){print$0}}' \
	 | awk '$3 = "exon"' | awk -v var="${start}" '$4 >=var' \
	 | awk -v var="${stop}" '$5 <=var'| awk '{print$12}' | uniq \
	 | sed 's/[^a-z  A-Z 0-9 .]//g' > transcripts.tmp ;
	echo "created transcripts.tmp"
	cat transcripts.tmp | while read tr_id; do
	echo $tr_id;
	echo ${PWD}/stringtie/$sample/$sample.gtf;
	#grep $tr_id tmp.gtf;
	grep -w $tr_id ${PWD}/stringtie/$sample/$sample.gtf >> ${PWD}/results/${gene}_$sample.gtf;
	done;
done;

############################################################################