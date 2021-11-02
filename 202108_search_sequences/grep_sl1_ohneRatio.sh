#!/bin/sh
set -e -x

SEQUENCES="sequences_comma.txt"

###### test 
# for i in /home/udo/worms/*.bam; do
	# echo ${i##*/} >> lat_1_test.txt
	# sequence="CCCAAGTTTGAGGAAACTCTTACACTCGACGAACAAAATAGTATGCACTGATGCGACGTAACAAAACGACTTATTCGTT"
	# counts=$(samtools view $i II | grep "$sequence" | wc -l)
	# echo "$counts" >> lat_1_test.txt
# done

########## count sequence in bam 
while IFS="," read -r field1 field2; do
	mkdir -p grep_output_ohneratio/
	if [ -f "grep_output_ohneratio/${field1}_readCount.txt" ]; then
			echo "grep_output_ohneratio/${field1}_readCount.txt exists."
		else
			touch "grep_output_ohneratio/${field1}_readCount.txt"
		fi
		echo "file \t count \t ratio" >> grep_output_ohneratio/${field1}_readCount.txt
	for i in /home/udo/worms/*.bam; do
		echo $i
		echo $field2
		string=$(echo "$field2" | tr -d ' ' | tr -d '[:blank:]' | tr -d '[:space:]') #remove all white spaces, blanks and spaces 
		counts=$(samtools view $i II | grep "$string" | wc -l) #count seq in mapping
		#total_counts=$(samtools view -c $i) #count reads 
		#ratio=echo "scale=5; $counts/$total_counts" | bc 	# calculate ratio
		echo "${i##*/}\t $counts" >> grep_output_ohneratio/${field1}_readCount.txt
	done
done < $SEQUENCES | awk -F',' '{print $1,$2}'
