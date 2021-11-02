#!/bin/sh
set -e -x

SEQUENCES="sequences_comma.txt"
BAMS="/home/udo/worms/"

###### test 
# for i in /home/biochemistry/2003KNO-0044/01.RawData/*.fastq.gz; do
	# seq_counts=$(zgrep -c 'ATGATGATG' $i ) #count appearance sequence in fastq 
	# total_counts=$(echo $(zcat $i | wc -l) / 4 | bc) #count all reads in fastq  
	# ratio_sequences=$(echo "scale=5;$seq_counts/$total_counts" | bc)
	# echo "${i##*/}\t $ratio_sequences" >> grep_from_fastq_test.txt
# done

###### WORKED

######## count reads containing sequences_1 test 
# for i in /home/biochemistry/2003KNO-0044/01.RawData/*.fastq.gz; do
	# seq_counts=$(zgrep -c 'ATGATGATG' $i ) #count appearance sequence in fastq 
	# total_counts=$(echo $(zcat $i | wc -l) / 4 | bc) #count all reads in fastq  
	# ratio_sequences=$(echo "scale=5;$seq_counts/$total_counts" | bc)
	# echo "${i##*/}\t $ratio_sequences" >> grep_from_fastq_test.txt
# done
###### WORKED

while IFS="," read -r field1 field2; do
	mkdir -p grep_output_fastq/
	echo $field2
	echo $field1
	touch grep_output_fastq/${field1}_grep_from_fastq.txt || exit
	for i in /home/biochemistry/2003KNO-0044/01.RawData/*.fastq.gz; do
		seq_counts=$(zgrep -c $field2 $i ) #count appearance sequence in fastq
		total_counts=$(echo $(zcat $i | wc -l) / 4 | bc) #count all reads in fastq
		ratio_sequences=$(echo "scale=5;$seq_counts/$total_counts" | bc)
		echo "${i##*/}\t $ratio_sequences" >> grep_output_fastq/${field1}_grep_from_fastq.txt 
	done
done < $SEQUENCES | awk -F',' '{print $1,$2}'