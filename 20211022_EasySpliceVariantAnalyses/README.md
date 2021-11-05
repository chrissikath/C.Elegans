# README.md (mapping, stringtie, lat-1, splice-variant-analysis)

Um die Analyse zu starten, lade den Ordner "20211022_EasySpliceVariantAnalyses" herunter, wo du die Analyse ausführen möchtest. 

## Part 1: splice_variants_analysis.sh
	
1.1 STAR INDEX </br>
1.2 Mapping  </br>
1.3 Sort bed files  </br>
2.0 Stingtie </br>
3.0 Extract Gene from gtf results </br>

Schritte:
- Programme installieren:
	- STAR version 2.7.4a
	- Stringtie version 1.3.3b
	- Beispielweise in einem environment installieren (z.B. conda)
- In der Datei splice-variant-analysis.sh Pfade/Gen-Angaben ändern, dann ausführen -> chmod +x splice-variant-analysis.sh -> ./splice-variant-analysis.sh </br>
	GENOME= .fa  </br>
	ANNO_GTF= .gtf </br>
	STARINDEX_OUTPUT = directory wo STARIndex hinsoll </br>
	RAW_SEQ= Ordner wo paired-end reads liegen (zippend) </br>
	SAMPLE_LIST= .txt -> welche samples genutzt werden (filename /t Endbezeichnung) </br>
	
	gene="Lat-1" # define gene name </br>
	start="8896841" # define start </br>
	stop="8908666"   # define stop </br>
	chr="II" #define chromosome </br>
	strand="+" #define strand </br>
		
	Example:
	nur der Name ohne _1 & _2, da paired-end sequencing
	``` G1-N2	d1_adult_gonad_wildtype_1
	G2-N2	d1_adult_gonad_wildtype_2
	G3-N2	d1_adult_gonad_wildtype_3
	G4-N2	d1_adult_gonad_wildtype_4
	G5-N2	d1_adult_gonad_wildtype_5 
	```
			
## Part 2: visualization_splice_variants.sh

> NOCH DRAN ALLE PACKAGES ZUM LAUFEN ZU BRINGEN 
> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1.0 Shell part (Vorbereitung fürs R script)  </br>
2.0 R Part (Visualisierung) </br>
	
- Programme installieren:
	- R version: 3.6.1 (conda create --name XYZ )
	```R
	conda install -c conda-forge r-base==3.6.1
	```
	- BiocManager
	
	```R
	if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")
	``` 
	- R-Packages installieren die in packages_functions_Alex.R benutzt werden </br>
		(DESeq2, edgeR, RColorBrewer, gplots , psych, VennDiagram, grDevices,
		calibrate, ggplot2, pheatmap, reshape, seqinr, Biostrings, plyr, taRifx)
	```R
	BiocManager::install(c("DESeq2", "edgeR", "RColorBrewer","gplots" , "psych", "VennDiagram"))
	BiocManager::install(c("grDevices", "calibrate", "ggplot2", "pheatmap", "reshape", "seqinr", "Biostrings", "plyr", "taRifx")) 
	```
		
- In der Datei Pfade .visualization_splice_variants.sh ändern, dann ausführen -> chmod +x visualization_splice_variants.sh -> ./visualization_splice_variants.sh
	- 1.0 files  </br>
		GENOME= .fa </br>
		packages_functions_script= packages_functions_Alex.R </br>
		exon_information_script= exon-information.r </br>
	- 2.0 files  </br>
		motifs_gene_file= motifs_gene  </br>
		set_working_dir= wo "results" hinsoll und alles liegt </br>
		visualization_splice_variants_script= visualization_splice_variants.R </br>
