#!/bin/sh
set -e -x

SEQUENCES="sequences_1.txt"
BAMS="/home/udo/worms/"

###### test 
#for i in /home/udo/worms/*.bam; do
#	echo ${i##*/} >> lat_1_test.txt
#	samtools view $i | grep "CTCAATTCATGCTGTCCTAGGAAACTATGGAAGATTCTCTGTCGCTGTTTG" | wc -l  >> lat_1_test.txt
#done
###### WORKED

#while IFS="," read -r field1 field2; do
#	mkdir -p grep_output/
#	echo $field2
#	touch grep_output/${field1}_readCount.txt || exit
#	for i in /home/udo/worms/*.bam; do
#		echo $i
##		echo ${i##*/} >> grep_output/${field1}_readCount.txt
#		samtools view $i "II" | grep $field2 | wc -l >> grep_output/${field1}_readCount.txt
#	done
#done < $SEQUENCES | awk -F',' '{print $1,$2}'

######## count reads containing sequences_1 test 
# for i in /home/biochemistry/2003KNO-0044/01.RawData/*.fastq.gz; do
	# seq_counts=$(zgrep -c 'ATGATGATG' $i ) #count appearance sequence in fastq 
	# total_counts=$(echo $(zcat $i | wc -l) / 4 | bc) #count all reads in fastq  
	# ratio_sequences=$(echo "scale=5;$seq_counts/$total_counts" | bc)
	# echo "${i##*/}\t $ratio_sequences" >> grep_from_fastq_test.txt
# done
## WORKED

while IFS="," read -r field1 field2; do
	mkdir -p grep_output_fastq/
	echo $field2
	echo $field1
	for i in /home/biochemistry/2003KNO-0044/01.RawData/*.fastq.gz; do
		if [ -f "grep_output_fastq/${field1}_readCount.txt" ]; then
			echo "grep_output_fastq/${field1}_readCount.txt exists."
		else
			touch "grep_output_fastq/${field1}_readCount.txt"
		fi
		total_counts=$(echo $(zcat $i | wc -l) / 4 | bc) #count all reads in fastq  
		seq_counts=$(zgrep -c ${field2} $i | cat ) #count appearance sequence in fastq 
		
		if [ "$seq_counts" -eq "0" ]; then
			echo "$seq_counts is 0 or some other arithmetic error occurred"
			ratio_sequences="NA"
		else
			ratio_sequences=$(echo "scale=5;$seq_counts/$total_counts" | bc)
		fi
		echo "${i##*/}\t $ratio_sequences" >> grep_output_fastq/${field1}_readCount.txt
	done
done < $SEQUENCES | awk -F',' '{print $1,$2}'
