#!/bin/sh
set -e -x

SEQUENCES="../data/sequences_comma.txt"
NOW=$( date '+%F' )
OUTPUT="./../analysis/search_sequences/${NOW}grep_output/"

###### test 
# for i in /home/udo/worms/*.bam; do
	# echo ${i##*/} >> lat_1_test.txt
	# sequence="CCCAAGTTTGAGGAAACTCTTACACTCGACGAACAAAATAGTATGCACTGATGCGACGTAACAAAACGACTTATTCGTT"
	# counts=$(samtools view $i II | grep "$sequence" | wc -l)
	# echo "$counts" >> lat_1_test.txt
# done

mkdir -p $OUTPUT #create output
while IFS="," read -r field1 field2; do
	if [ -f "${OUTPUT}${field1}_readCount.txt" ]; then
			echo "${OUTPUT}${field1}_readCount.txt exists."
		else
			touch "${OUTPUT}${field1}_readCount.txt"
		fi
		#echo "file \t count \t ratio" >> ${OUTPUT}${field1}_readCount.txt
	for i in /home/udo/worms/*.bam; do
		echo $i
		echo $field2
		string=$(echo "$field2" | tr -d ' ' | tr -d '[:blank:]' | tr -d '[:space:]') #remove all white spaces, blanks and spaces 
		counts=$(samtools view $i II | grep "$string" | wc -l) #count seq in mapping
		total_counts=$(samtools view -c $i) #count reads 
		if [ "$seq_counts" -eq "0" ]; then
			echo "$seq_counts is 0 or some other arithmetic error occurred"
			ratio="0"
		else
			ratio=$(echo "$counts/$total_counts" | bc -l)
		fi
		echo "${i##*/}\t $counts \t $ratio" >> ${OUTPUT}${field1}_readCount.txt
	done
done < $SEQUENCES | awk -F',' '{print $1,$2}'
