# README.md (mapping, stringtie, lat-1, splice-variant-analysis)

## Part 1: splice_variants_analysis.sh
	
1.1 STAR INDEX <tr>
1.2 Mapping
1.3 Sort bed files
2.0 Stingtie
3.0 Extract Gene from gtf results
	
- Versionen:
	- STAR version 2.7.4a
	- Stringtie version 1.3.3b
- In der Datei Pfade ändern, dann ausführen -> ./splice-variant-analysis.sh 
	- Extra-Dateien:
		GENOME= .fa 
		ANNO_GTF= .gtf
		RAW_SEQ= Ordner wo paired-end reads liegen (zippend)
		SAMPLE_LIST= .txt -> welche samples genutzt werden (filename /t Endbezeichnung)
		Example: 	G1-N2	d1_adult_gonad_wildtype_1
					G2-N2	d1_adult_gonad_wildtype_2
					G3-N2	d1_adult_gonad_wildtype_3
					G4-N2	d1_adult_gonad_wildtype_4
					G5-N2	d1_adult_gonad_wildtype_5
			
## Part 2: visualization_splice_variants.sh

1.0 Shell part (Vorbereitung fürs R script)
2.0 R Part (Visualisierung)
	
- Versionen: 
	- R version: 3.4.2
	- BiocManager::install(version = "3.10")
		R-Packages installieren die in packages_functions_Alex.R benutzt werden
		(DESeq2, edgeR, RColorBrewer,gplots , psych,VennDiagram,grDevices,
		calibrate, ggplot2, pheatmap, reshape, seqinr, Biostrings, plyr, taRifx)
	- In der Datei Pfade ändern, dann ausführen -> ./visualization_splice_variants.sh
	- #1.0 files
	GENOME= .fa
	R_scripts= packages_functions_Alex.R
	exon_information_script= exon-information.r
		- #2.0 files
	motifs_adgrl1_file= motifs_adgrl1 
	packages_functions_script= packages_functions_Alex.R
	set_working_dir= wo "resuls" hinsoll und alles liegt
	visualization_splice_variants_script= visualization_splice_variants.R
