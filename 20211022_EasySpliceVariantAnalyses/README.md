# README.md (mapping, stringtie, lat-1, splice-variant-analysis)

## Part 1: splice_variants_analysis.sh
	
1.1 STAR INDEX </br>
1.2 Mapping  </br>
1.3 Sort bed files  </br>
2.0 Stingtie </br>
3.0 Extract Gene from gtf results </br>
	
- Versionen:
	- STAR version 2.7.4a
	- Stringtie version 1.3.3b
- In der Datei Pfade ändern, dann ausführen -> ./splice-variant-analysis.sh 
	- Extra-Dateien: </br>
		GENOME= .fa  </br>
		ANNO_GTF= .gtf </br>
		RAW_SEQ= Ordner wo paired-end reads liegen (zippend) </br>
		SAMPLE_LIST= .txt -> welche samples genutzt werden (filename /t Endbezeichnung) </br>
		Example: 	
		''' 	G1-N2	d1_adult_gonad_wildtype_1
			G2-N2	d1_adult_gonad_wildtype_2
			G3-N2	d1_adult_gonad_wildtype_3
			G4-N2	d1_adult_gonad_wildtype_4
			G5-N2	d1_adult_gonad_wildtype_5 '''
			
## Part 2: visualization_splice_variants.sh

1.0 Shell part (Vorbereitung fürs R script)  </br>
2.0 R Part (Visualisierung) </br>
	
- Versionen: 
	- R version: 3.4.2
	- BiocManager::install(version = "3.10") </br>
		R-Packages installieren die in packages_functions_Alex.R benutzt werden </br>
		(DESeq2, edgeR, RColorBrewer,gplots , psych,VennDiagram,grDevices,
		calibrate, ggplot2, pheatmap, reshape, seqinr, Biostrings, plyr, taRifx)
	- In der Datei Pfade ändern, dann ausführen -> ./visualization_splice_variants.sh
	- #1.0 files  </br>
		GENOME= .fa </br>
		R_scripts= packages_functions_Alex.R </br>
		exon_information_script= exon-information.r </br>
	- #2.0 files  </br>
		motifs_adgrl1_file= motifs_adgrl1  </br>
		packages_functions_script= packages_functions_Alex.R </br>
		set_working_dir= wo "resuls" hinsoll und alles liegt </br>
		visualization_splice_variants_script= visualization_splice_variants.R </br>
