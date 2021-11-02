# R script for visualization
args <- commandArgs()
print(args)

if (length(args)==0){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==1){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==2){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==3){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==4){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==5){
	stop("You need 4 args: gene_name, motifs_adgrl1_file, packages_functions_script ,set_working_dir ")
}else if (length(args)==6){
	gene.name <- args[3]	
	motifs_adgrl1_file <- args[4]
	packages_functions_script <- args[5]
	set_working_dir <- args[6]
}

source(file = packages_functions_script)
setwd (dir = set_working_dir)

# import motif table

motif_table <- read.table(file = motifs_adgrl1_file, sep = "\t", stringsAsFactors = FALSE)

# build the transcripts, analyse sequence, find and analyse ORF

transcripts <- GetORFAndProteinSeqFromConcatGTFAndMotifFinder(
    motif_table = motif_table,
    exons_fasta = "exon_seqs.fa", 
    conc_tr.gtf = "all_transcripts.gtf", 
    pdfname = "sequence_and_orf_analysis.pdf", 
    orfseqname = "orfs.fa", 
    protseqname = "peptides.fa")

#dev.off() # just in case
write.table(x = transcripts, file = paste("transcripts_", gene.name, '.txt', sep = ''), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

transcripts <- read.table(file = paste("transcripts_", gene.name, '.txt', sep = ''), 
                          stringsAsFactors = FALSE, sep = "\t", header = FALSE)
						  # transcript plot

exons <- read.table(file = "results/unique_exons.gtf", sep = "\t", stringsAsFactors = FALSE)

domains <- c("S-Pep","GPS","TM1")

domain_locations <- c(
    paste( 7.26,7.81,sep=" "),
    paste(20.79,20.90,sep=" "),
    paste(20.93,20.98,sep=" ")
    )
    
domain_df <- as.data.frame(cbind(domains,domain_locations))

domain_df <- as.data.frame(cbind(domains,domain_locations))
domain_df <- FALSE # dont print motifs at transcripts.pdf

pdf("Transcript_variants.pdf", width = 18, height = 10)

PlotAllTranscriptsInExonGraph(
    exon_data = exons,
    transcripts = transcripts,
    #transcript_rows = c(1:9),
    gene_name = gene.name,
    vis.ref.file = "vis.ref.file.l1",
    domains = domain_df)
	
dev.off()

#################################DONE#################################################
