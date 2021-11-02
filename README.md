# C elegans analysis for Victoria and Daniel
Splice Variant Analysis for LAT-1

## Victoria 

### Exon quantification

- got exon positions from same genome annotation like mappings 
	../data/chromosom2_exons.bed
- count reads hitting exons with bathometer 
	../scripts/bathometer.sh
- *output*: FPKM values for each exon 
	../analysis/exon_quantification/exon_quantification.tsv
- File hochgeladen in webtool
- Exon-Name include Chromosom:Start-Stop , Strand, ParentTranscript
- "Parent-Transcipt" -> hierarchical because one CDS/transcript (parent) can have more exons (child) 

### Search specific sequences (found in RACE)

- count for each sequence and its reversed sequence 
	../data/sequence_and_reverse.txt
- normalise for read-coverage, proportion of these reads per sample:
	- ratio = count/allReads

## Daniel 	

### genes bathometer 

- run bathometer for given bed file with genes 
	../analysis/20210816_genes_bathometer/run_bathometer.sh
- file hochgeladen in webtool 


### Figure LAT-1 splice variantes 

- checked reproducability of Alex results: 

#### worms

1.0 splice_variants_analysis.sh 
- Step 1.1: StarIndex 
- Step 1.2: Mapping with Star -> bam
- Step 1.2: Sort -> sorted.bam
- Step 2.0: Quantification with Stringtie -> gtf 
- Step 3.0. Extract wanted gene (LAT-1) from results -> lat1.gtf

2.0 visualization_splice_variants.sh
- runs R script -> PFADE ÄNDERN


### Fastqs per transcript result for cloning
- gtf2bed -> .bed -> sort for gene_id > getfasta from .fa > .fasta 
``` gtf2bed < lat1_d1_adult_gonad_wildtype_2.gtf > lat1_d1_adult_gonad_wildtype_2.bed  ```
```cat lat1_d1_adult_whole_wildtype_4.bed | sort -k10 > lat1_d1_adult_whole_wildtype_4_sorted.bed ```
``` bedtools getfasta -fi ../../../data/wbps14/caenorhabditis_elegans.PRJNA13758.WBPS14.genomic.fa -bed lat1_d1_adult_gonad_wildtype_5_sorted.bed -fo lat1_d1_adult_gonad_wildtype_5.fasta -name+ ```