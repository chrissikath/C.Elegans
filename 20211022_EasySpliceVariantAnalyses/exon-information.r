### R script for exon information
### 
### script will take input of 
### 1, gene name 
### 2, sample names 
### assuming that input files are stored in the current working directory


if (interactive() == FALSE) {
	Args <- commandArgs()
	if (length(Args) == 2)  {
		stop('please provide the sample names to be loaded')
	}
}

# R
# Args <- c('R', 'vanilla', 'radish', 'adult_male_1', 'adult_male_2')

source(file = "/home/christina/C_elegans/analysis/20210906_AlexRedone/scripts/packages_functions_Alex.R")
getwd()

# print gene of interest to provide possibillity to interupt
gene.name <- Args[3]
print(paste('gene of interest:', gene.name))

# create an empty list for all gtf files
gtf.list <- list()
merged.exons <- FALSE

for (arg in Args[4:length(Args)]) {
	if (grepl('.gtf', arg)) {
		sample.name <- strsplit(arg, '\\.')[[1]][[1]]
	}
	else {
		sample.name <- arg
	}
	print(paste('sample: ', sample.name, sep = ''))
	# for each sample: load gtf file, append list with it and name it according to sample
	new.file <- read.table(file = paste(sample.name, '.gtf', sep = ''), sep = '\t', stringsAsFactors = FALSE)
	gtf.list <- append(gtf.list, list(new.file))

	names(gtf.list)[length(gtf.list)] <- sample.name
	
	# get unique exons for each file
	uniq.exon <- MergeDuplicatesInGtfWithFreq(new.file)
	if (typeof(merged.exons) == 'logical') {
		merged.exons <- uniq.exon
		# retrieve data for chromosome and strand from gtf file 
		# also check, whether all entries in gtf come from same chromosome and strand 
		chromosome <- unique(new.file[, 1])
		strand <- unique(new.file[, 7])
		if (length(chromosome) > 1 | length(strand) > 1) {
			print('WARNING: Check location - more than 1 strand and/or chromosome')
		} 
	}
	else {
		merged.exons <- rbind(merged.exons, uniq.exon)
		if (length(unique(new.file[, 1])) > 1 || unique(new.file[, 1]) != chromosome) { 
			print('WARNING: Check location - more than 1 chromosome')
		}
		if (length(unique(new.file[, 7])) > 1 || unique(new.file[, 7]) != strand ) {
			print('WARNING: Check location - more than 1 strand')
		} 

	}
}

print(str(gtf.list))

merged.exons <- MergeGtfInformationWithFreq(merged.exons,(length(Args) - 2))
ExportExonDataFrameAsGtf(merged.exons, chromosome, strand, gene.name, "unique_exons.gtf")

for (sample in 1:length(gtf.list)) {
	AppendCovAndUniqueIdentifierAndExport(gtf = gtf.list[[sample]], meta_data = merged.exons, 
		filename = paste(names(gtf.list)[sample], '_append.gtf', sep = ''))
}
