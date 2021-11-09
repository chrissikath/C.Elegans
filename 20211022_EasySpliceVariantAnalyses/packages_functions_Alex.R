###########################################################################################
###                                                                                     ###
###             R Script for easy loading of packages and often used functions          ###
###                                                                                     ###
###########################################################################################

###########################################################################################
###                                 0. How to use!                                      ###
###########################################################################################

# Implement following code at the beginning of your R Script.
# source(file="/home/user/mpi/mnt/raid01/data/Alex/scripts/packages_functions_Alex.R")

###########################################################################################
###                                1. Load packages!                                    ###
###########################################################################################

#source("https://bioconductor.org/biocLite.R") # for easy installing of BioConductor packages

# BiocManager::install(c("DESeq2", "edgeR, RColorBrewer,gplots , psych, GMD,VennDiagram,grDevices,\
# venneuler,   calibrate, ggplot2, pheatmap, reshape, heatmap.plus, Heatplus, seqinr, Biostrings, plyr, taRifx))
#library(DESeq2) # statistical analysis for RNA-Seq
#library(edgeR)  # statistical analysis for RNA-Seq
#library(RColorBrewer) # for color palets
#library(gplots) # plotting
#library(psych) # stats
# library(GMD) # stats
#library(VennDiagram)
#library(grDevices)		# needed to save the venn diagrams as pdf
# library(venneuler)
#library(calibrate) # plotting
#library(ggplot2) # plotting
#library(pheatmap) # = pretty heatmap - plotting
#library(reshape) # for some basic functions
# library(heatmap.plus) # plotting
# library(Heatplus) # plotting
library(seqinr) # for open reading frame detection & genetic code
#library(Biostrings) # for open reading frame detection & genetic code
#library(plyr) # for rounding functions
#library(taRifx) # for remove.factors() from a dataframe function
# install.packages("")

###########################################################################################
###                                2. Load functions!                                   ###
###########################################################################################


### Table of content ###

# I General functions #

#   - "zscore"
#     Normalisation function. Input: String of numeric arguments. Calculates values minus mean & divides by standard deviation.
#   - "TestAColorPalette"
#     Visual function. Opens x11 window with a pie chart where each pie is colored according to color palette input.
#   - "alex11"
#     Open x11 - window with default size which fits to bioscho03 main screen. (Pun intended)
#   - "heatmap.3.custom"
#     Visual function. Plot a heatmap. (Function was found by vera online.)
#   - "uncollapse"
#     Uncollapse a string. It is advised to append [[1]] directly after this command since it returns a list.
#   - "allrows" + "allcols"
#     Get vector with all positive integers from 1 to the length of the rows of a dataframe. (Good for loop)
#   - "RemoveElements"
#     Insert character and remove number of characters (from the start or from the end)
#   - "ReadIn"
#     Shorter read.table with sep = "\t" and all specifications I usually have

# II NGS-specific functions # 

#   - "CalculatePeptideMass"
#     Calculate the mass of all residues in a protein sequence.
#   - "findPotentialStartsAndStops"
#     Finds start and stop codons. (Origin: http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter7.html)
#   - "findPotentialStartsAndStops2"
#     Faster version of "findPotentialStartsAndStops". (Credits: http://www.talkstats.com/showthread.php/28335-Finding-start-and-stop-codons-in-a-DNA-sequence)
#   - "plotPotentialStartsAndStops"
#     Plots the output of findPotentialStartsAndStops.
#   - "plotPotentialStartsAndStopsChangeX"
#     Same function as plotPotentialStartsAndStops, but x-axis is changeable (useful for comparisons of multiple sequences of different length).
#   - "findORFsinSeq"
#     Finds ORFs in a nucleotide sequence.  
#   - "plotORFsinSeq"
#     Plot the output of findORFsinSeq.
#   - "plotORFsinSeqChangeX"
#     Same function as plotORFsinSeq, but x-axis is changeable (useful for comparisons of multiple sequences of different length).
#   - "RemoveLowlyAbundantTranscripts"
#     Apply exon or transcript wise cutoff to loaded gtf file.
#   - "LoadGTFAndRemoveLowlyAbundantTranscripts"
#     Load GTF file and apply exon or transcript cutoff.
#   - "MergeDuplicatesInGtf"
#     Extract a non-redundant exon dataset including the coverage of a gtf file. 
#   - "MergeDuplicatesInGtfWithFreq"
#     Same function as "MergeDuplicatesInGtf", but with added frequency of each exon.
#   - "MergeGtfInformation"
#     Merge multiple samples preprocessed by "MergeDuplicatesInGtf".
#   - "MergeGtfInformationWithFreq"
#     Same as "MergeGtfInformation", but with added frequency of each exon.
#   - "MergeGtfInformationWithFreqInclCondition"
#     Same as "MergeGtfInformationWithFreq", but with condition-wise coverage values.
#   - "ExportDataFrameAsGtf"
#     Export output from "MergeGtfInformation*" to IGV-compatible GTF-file.
#   - "ExportDataFrameAsGTFINclCondition"
#     Same as "ExportDataFrameAsGtf", but includes condition-wise coverage values (made from "MergeGtfinfromationWithFreqInclCondition").
#   - "AppendCovAndUniqueIdentifierAndExport"
#     Add meta-information from multiple merged samples about exons (coverage & unique identifier) to a transcript gtf file.
#   - "ConvertGTFtoGFFandExport"
#     Convert transcript gtf file to an gff file based on meta.data about exons. Was needed for DEXSeq.
#   - "GetORFAndProteinSeqFromConcatGTF"
#     Assemble transcript sequences from exons.fa & concat transcript.gtf, find & plot ORFs, translate & output.
#   - "plotMotifsInSequence"
#     Finds (exact) amino acid motifs/sequences in a protein sequence and plot it.
#   - "GetORFAndProteinSeqFromConcatGTFAndMotifFinder"
#     "GetORFAndProteinSeqFromConcatGTF" + motif finder + peptide mass calculator.
#   - "MergeDuplicateORFs"
#     Find same ORFs and group the transcripts with same ORFs.
#   - "MergeDuplicateTranscripts"
#     Find same transcripts (based on exon composition). Define start- and endgroup when having different 5' / 3' - UTRs.
#   - "PlotTranscriptsInExonGraph"
#     Exon plot similar to DEXSeq with dynamic genetic region and shortened introns.
#   - "PlotTranscriptsInExonGraphInclCondition"
#     Same function as "PlotTranscriptsInExonGraph" with color coding based on condition coverage.
#   - "GetMatsResults"
#     Load rMATS results to R workspace.
#   - "RemoveMatsResFromEnv"
#     Remove rMATS results from R workspace. STRONGLY recommended when using "GetMatsResults" repetitively for different genes.
#   - "CleanUpStringtieOutput"
#     Outdated function! Was used for automatisation of the first steps of the Stringtie+R pipeline. Does not work for conditional coverage values!
#   - "PlotMatsInExonGraph*"
#     Several functions for combining the rMATS output and the "PlotTranscriptInExonGraph*" functions. Scroll to functions for further information.
#   - "GetRefTranscripts"
#     Give GTF file with reference transcripts as input and overwrite the exons with your unique ID.
#   - "Translate3Frames"
#     Input: Nucleotide sequence | Output: Peptide Sequence table (3 frames)
#   - "PlotExonsWithoutOverlap"
#     Input: vis.ref.file from PlotTrinExonGraph | Output: X-Axis of PlotTranscriptInExonGraph and below the compositions of this
#   -"PlotExonCountDistribution"
#     Input: vis.ref.file & bedtools-multicov-output-file | Output: Plot with exon coverage
#   -"TripleLength"
#     Input: Character (Example: Peptide sequence) | Output: 3 times the amount of parts of the character (the number of bases for the peptide sequence)
#   -"DoPeptideAndNucleotideSequenceMatch"
#     Input: Peptide Sequence and DNA-Sequence | Output: Test result if DNA-Sequence is coding Peptide Sequence
#   -"ReverseComplement"
#     Own variation of BioStrings "reverseComplement" function without annoying BioString object formats
#   -"ImportReferenceFasta"
#     Give a file which is a fasta, which means, that it displays a long sequence in line-blocks of 50 characters. Output: A table with the concatenated sequence string and the self-given name.
#   -"AnalyseTranscripts"
#     Analyse tr_-table regarding frequency in samples and groups (tissues)
#   -"AddSamplesAndRelFPKMToTranscripts"
#     Split function "AnalyseTranscripts" I: Analyse the samples and groups from tr_d1 only and calculate the fpkm_wise percentage of each transcript
#   -"CompareTranscriptsOverGroups"
#     Split function "AnalyseTranscripts" II: Analyse the output of the function above to find out how often each ORF exists in which sample / tissue + calculate it
#   -"CompareReferencesAndORFs"
#     Input: Output of function above and Output of ImportReferenceFasta
#   -"VerboseInterpretation"
#     Input: Output of function above 
#   -"ComparativeTranscriptPlot"
#     Similar to PlotTrInExonGraph, but with ComparedTranscripts and a Heatmap with tissue/type quantification
#   -"ExonCoveragePlot"
#     Input: Samtools Depth Output file --> Output is a plot which shows the exonic coverage
#   -"ReadGTFApplyCutoff"
#     Load GTF and apply transcript-based cutoff: supported options: relative or absolute | FPKM, TPM & cov
#   -"CutOffORFs"
#     Based on fraction in every tissue, a cutoff is being made.

########################

# I General functions # 

TestAColorPalette <- function(colors){
    cur <- dev.cur() # get active device
    x11() # open new
    pie(rep(1, length(colors)), labels = sprintf("%d (%s)", seq_along(colors), colors), col = colors) # plot pie with all colors
    invisible(dev.set(cur)) # change back to active device
}

alex11 <- function(width = 19, heigth = 10){
    x11(width = width, height = heigth)
}

heatmap.3.custom <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, distfun = dist, hclustfun = hclust, dendrogram = c("both","row", "column", "none"), symm = FALSE,
                             scale = c("none","row", "column"), na.rm = TRUE, revC = identical(Colv,"Rowv"), add.expr, breaks, symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                             col = "heat.colors", colsep, rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("none", "column","row", "both"),
                             tracecol = "cyan", hline = median(breaks), vline = median(breaks), linecol = tracecol, margins = c(5,5), ColSideColors, RowSideColors,
                             side.height.fraction=0.3, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, key = TRUE,keysize = 1.5, density.info = c("none", "histogram", "density"),
                             denscol = tracecol, symkey = max(x < 0, na.rm = TRUE) || symbreaks,densadj = 0.25, main = NULL,
                             xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL,lwid = NULL, NumColSideColors = 1,NumRowSideColors = 1,KeyValueName="Value",...){
    
    invalid <- function (x) {
        if (missing(x) || is.null(x) || length(x) == 0)
            return(TRUE)
        if (is.list(x))
            return(all(sapply(x, invalid)))
        else if (is.vector(x))
            return(all(is.na(x)))
        else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                     c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                     c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                          length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        
        if (!missing(ColSideColors)) {
            #if (!is.matrix(ColSideColors))
            #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
            lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }
        
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(colnames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
    
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
              col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                       lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
             col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
              xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                  lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                  col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

uncollapse <- function(string){
    strsplit(as.character(string), split = " ")
}

allrows <- function(dataframe){
    rows <- c(1:length(rownames(dataframe)))
    return(rows)
}

allcols <- function(dataframe){
    cols <- c(1:length(colnames(dataframe)))
    return(cols)
}

RemoveElements <- function(char, amount = 1, start = "last"){
    if(start == "last"){
        char <- as.character(char)
        out <- substr(char,1, nchar(char)-amount)
    }else if(start == "first"){
        len <- length(s2c(char))
        min <- amount + 1
        out <- c2s(s2c(char)[min:len])
        
    }else{print("Pleasy specify where to start - default: last, other option: first")}
    return(out)
    
}

ReadIn <- function(file){
    file <- read.table(file = file, stringsAsFactors = FALSE, header = FALSE, sep = "\t")
    return(file)
}


############################

# II NGS-specific funtions # 

z_score <- function(x){
    
    mean <- sum(x)/length(x)
    sd <- sd(x)
    result <- (x-mean)/sd
    return(result)
}

CalculatePeptideMass <- function(peptide_sequence, unit){ 
    # dont forget H @ N & OH @ C  Terminus
    library(seqinr)
    # define residue mass (residue = amino acid in peptide chain)
    # source: http://www.matrixscience.com/help/aa_help.html [Juli,13th,2017] - took the average mass
    a <- 71.0779
    r <- 156.1857
    n <- 114.1026
    d <- 115.0874
    c <- 103.1439
    e <- 129.114
    q <- 128.1292
    g <- 57.0513
    h <- 137.1393
    i <- 113.1576
    l <- 113.1576
    k <- 128.1723
    m <- 131.1961
    f <- 147.1739
    p <- 97.1152
    s <- 87.0773
    t <- 101.1039
    u <- 150.0379
    w <- 186.2099
    y <- 163.1733
    v <- 99.1311
    
    mass_index <- data.frame(c(a,r,n,d,c,e,q,g,h,i,l,k,m,f,p,s,t,u,w,y,v),row.names = c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","U","W","Y","V"))
    colnames(mass_index) <- c("mass")
    
    # sequence must be upper case, since mass value is stored in lower case
    
    prot_seq <- toupper(peptide_sequence)
    
    # split into single characters
    
    prot_seq <- s2c(prot_seq)
    
    mass <- 0
    for (i in prot_seq[-length(prot_seq)]){ #remove last *
        #print(i)
        mass <- mass + mass_index[i,"mass"]
        #print(mass)
    }
    nucl <- 0       # for safety: count how many possible Nucleotides (ATGC) could be counted
    for (i in prot_seq[-length(prot_seq)]){
        if(i == "C"){
            nucl <- nucl+1
        }else if(i == "A"){
            nucl <- nucl +1
        }else if (i == "T"){
            nucl <- nucl +1
        }else if (i == "G"){
            nucl <- nucl+1
        }
    }
    
    if ((length(prot_seq)-nucl)<=(length(prot_seq)*0.1)){ # if the "nucleotides" form more than 90% of the sequence, print an error message
        print("Careful: >90% of Sequence consists of A,T,G & C. Proceed with caution, propably a nucleotide sequence was used as input.")
    }
    
    if (unit == "Da" ){
        #print(paste("Mass of residues in Dalton:---  ",mass,"  ---.",sep = ""))
        return(mass)
    }
    if(unit == "kDa"){
        #print(paste("Mass of residues in kDa:---   ",mass/1000,"   ---.", sep = ""))
        return(mass/1000)
    }
    if(unit == "MDa"){
        #print(paste("Mass of residues in MDa:---   ",mass/1000000,"   ---.", sep =""))
        return(mass/1000000)
    }
    
}

findPotentialStartsAndStops <- function(sequence){
    # Define a vector with the sequences of potential start and stop codons
    codons            <-c("atg", "taa", "tag", "tga")
    # Find the number of occurrences of each type of potential start or stop codon
    for (i in 1:4)
    {
        codon <- codons[i]
        # Find all occurrences of codon "codon" in sequence "sequence"
        occurrences <- matchPattern(codon, sequence)
        # Find the start positions of all occurrences of "codon" in sequence "sequence"
        codonpositions <- start(occurrences)
        # Find the total number of potential start and stop codons in sequence "sequence"
        numoccurrences <- length(codonpositions)
        if (i == 1)
        {
            # Make a copy of vector "codonpositions" called "positions"
            positions <- codonpositions
            # Make a vector "types" containing "numoccurrences" copies of "codon"
            types <- rep(codon, numoccurrences)
        }
        else
        {
            # Add the vector "codonpositions" to the end of vector "positions":
            positions   <- append(positions, codonpositions, after=length(positions))
            # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
            types       <- append(types, rep(codon, numoccurrences), after=length(types))
        }
    }
    # Sort the vectors "positions" and "types" in order of position along the input sequence:
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
    # Return a list variable including vectors "positions" and "types":
    mylist <- list(positions,types)
    return(mylist)
}

findPotentialStartsAndStops2 <- function(sequence, codons = c("atg", "taa", "tag", "tga")){
    # find the starting position of the matches for all the codons
    startpositions <- sapply(codons, function(x){start(matchPattern(x, sequence))})
    # combine them into a vector and sort by position
    output <- sort(do.call(c, startpositions))
    
    # If we really *need* a list
    #output <- list(as.numeric(output), names(output))
    
    return(output)
    
}

plotPotentialStartsAndStops <- function(sequence){
    # Define a vector with the sequences of potential start and stop codons
    codons <- c("atg", "taa", "tag", "tga")
    # Find the number of occurrences of each type of potential start or stop codon
    for (i in 1:4){
        codon <- codons[i]
        # Find all occurrences of codon "codon" in sequence "sequence"
        occurrences <- matchPattern(codon, sequence)
        # Find the start positions of all occurrences of "codon" in sequence "sequence"
        codonpositions <- start(occurrences)
        # Find the total number of potential start and stop codons in sequence "sequence"
        numoccurrences <- length(codonpositions)
        if (i == 1)
        {
            # Make a copy of vector "codonpositions" called "positions"
            positions   <- codonpositions
            # Make a vector "types" containing "numoccurrences" copies of "codon"
            types       <- rep(codon, numoccurrences)
        }
        else
        {
            # Add the vector "codonpositions" to the end of vector "positions":
            positions   <- append(positions, codonpositions, after=length(positions))
            # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
            types       <- append(types, rep(codon, numoccurrences), after=length(types))
        }
    }
    # Sort the vectors "positions" and "types" in order of position along the input sequence:
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
    # Make a plot showing the positions of the start and stop codons in the input sequence:
    # Draw a line at y=0 from 1 to the length of the sequence:
    x  <- c(1,nchar(sequence))
    y <- c(0,0)
    plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide # in sequence", ylab="Reading frame",
         main="Predicted start (green) and stop (red) codons")
    segments(1,1,nchar(sequence),1)
    segments(1,2,nchar(sequence),2)
    # Add the x-axis at y=0:
    axis(1, pos=0)
    # Add the y-axis labels:
    text(0.9,0.5,"+1")
    text(0.9,1.5,"+2")
    text(0.9,2.5,"+3")
    # Draw in each predicted start/stop codon:
    numcodons <- length(positions)
    for (i in 1:numcodons)
    {
        position <- positions[i]
        type <- types[i]
        remainder <- (position-1) %% 3
        if    (remainder == 0) # +1 reading frame
        {
            if (type == "atg") { segments(position,0,position,1,lwd=2,col="mediumseagreen") }
            else               { segments(position,0,position,1,lwd=1,col="red")}
        }
        else if (remainder == 1)
        {
            if (type == "atg") { segments(position,1,position,2,lwd=2,col="mediumseagreen") }
            else               { segments(position,1,position,2,lwd=1,col="red")}
        }
        else if (remainder == 2)
        {
            if (type == "atg") { segments(position,2,position,3,lwd=2,col="mediumseagreen") }
            else               { segments(position,2,position,3,lwd=1,col="red")}
        }
    }
}

plotPotentialStartsAndStopsChangeX <- function(sequence, xscale){
    # Define a vector with the sequences of potential start and stop codons
    codons <- c("atg", "taa", "tag", "tga")
    # Find the number of occurrences of each type of potential start or stop codon
    for (i in 1:4){
        codon <- codons[i]
        # Find all occurrences of codon "codon" in sequence "sequence"
        occurrences <- matchPattern(codon, sequence)
        # Find the start positions of all occurrences of "codon" in sequence "sequence"
        codonpositions <- start(occurrences)
        # Find the total number of potential start and stop codons in sequence "sequence"
        numoccurrences <- length(codonpositions)
        if (i == 1)
        {
            # Make a copy of vector "codonpositions" called "positions"
            positions   <- codonpositions
            # Make a vector "types" containing "numoccurrences" copies of "codon"
            types       <- rep(codon, numoccurrences)
        }
        else
        {
            # Add the vector "codonpositions" to the end of vector "positions":
            positions   <- append(positions, codonpositions, after=length(positions))
            # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
            types       <- append(types, rep(codon, numoccurrences), after=length(types))
        }
    }
    # Sort the vectors "positions" and "types" in order of position along the input sequence:
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
    # Make a plot showing the positions of the start and stop codons in the input sequence:
    # Draw a line at y=0 from 1 to the length of the sequence:
    x  <- c(1,xscale)           ## changed it to xscale (parameter given at the beginning)
    y <- c(0,0)
    plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide # in sequence", ylab="Reading frame",
         main="Predicted start (green) and stop (red) codons")
    segments(1,1,nchar(sequence),1)
    segments(1,2,nchar(sequence),2)
    # Add the x-axis at y=0:
    axis(1, pos=0)
    # Add the y-axis labels:
    text(0.9,0.5,"+1")
    text(0.9,1.5,"+2")
    text(0.9,2.5,"+3")
    # Draw in each predicted start/stop codon:
    numcodons <- length(positions)
    for (i in 1:numcodons)
    {
        position <- positions[i]
        type <- types[i]
        remainder <- (position-1) %% 3
        if    (remainder == 0) # +1 reading frame
        {
            if (type == "atg") { segments(position,0,position,1,lwd=2,col="mediumseagreen") }
            else               { segments(position,0,position,1,lwd=1,col="red")}
        }
        else if (remainder == 1)
        {
            if (type == "atg") { segments(position,1,position,2,lwd=2,col="mediumseagreen") }
            else               { segments(position,1,position,2,lwd=1,col="red")}
        }
        else if (remainder == 2)
        {
            if (type == "atg") { segments(position,2,position,3,lwd=2,col="mediumseagreen") }
            else               { segments(position,2,position,3,lwd=1,col="red")}
        }
    }
}

findORFsinSeq <- function(sequence){
    require(Biostrings)
    # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
    mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
    orfstarts <- numeric()
    orfstops <- numeric()
    # Make a vector "orflengths" to store the lengths of the ORFs
    orflengths <- numeric()
    # Print out the positions of ORFs in the sequence:
    # Find the length of vector "positions"
    numpositions <- length(positions)
    # There must be at least one start codon and one stop codon to have an ORF.
    if (numpositions >= 2)
    {
        for (i in 1:(numpositions-1))
        {
            posi <- positions[i]
            typei <- types[i]
            found <- 0
            while (found == 0)
            {
                for (j in (i+1):numpositions)
                {
                    posj  <- positions[j]
                    typej <- types[j]
                    posdiff <- posj - posi
                    posdiffmod3 <- posdiff %% 3
                    # Add in the length of the stop codon
                    orflength <- posj - posi + 3
                    if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
                    {
                        # Check if we have already used the stop codon at posj+2 in an ORF
                        numorfs <- length(orfstops)
                        usedstop <- -1
                        if (numorfs > 0)
                        {
                            for (k in 1:numorfs)
                            {
                                orfstopk <- orfstops[k]
                                if (orfstopk == (posj + 2)) { usedstop <- 1 }
                            }
                        }
                        if (usedstop == -1)
                        {
                            orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                            orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                            orflengths <- append(orflengths, orflength, after=length(orflengths))
                        }
                        found <- 1
                        break
                    }
                    if (j == numpositions) { found <- 1 }
                }
            }
        }
    }
    # Sort the final ORFs by start position:
    indices <- order(orfstarts)
    orfstarts <- orfstarts[indices]
    orfstops <- orfstops[indices]
    # Find the lengths of the ORFs that we have
    orflengths <- numeric()
    numorfs <- length(orfstarts)
    for (i in 1:numorfs)
    {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        orflength <- orfstop - orfstart + 1
        orflengths <- append(orflengths,orflength,after=length(orflengths))
    }
    mylist <- list(orfstarts, orfstops, orflengths)
    return(mylist)
}

plotORFsinSeq <- function(sequence){
    # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
    mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
    orfstarts <- numeric()
    orfstops <- numeric()
    # Make a vector "orflengths" to store the lengths of the ORFs
    orflengths <- numeric()
    # Print out the positions of ORFs in the sequence:
    numpositions <- length(positions) # Find the length of vector "positions"
    # There must be at least one start codon and one stop codon to have an ORF.
    if (numpositions >= 2)
    {
        for (i in 1:(numpositions-1))
        {
            posi <- positions[i]
            typei <- types[i]
            found <- 0
            while (found == 0)
            {
                for (j in (i+1):numpositions)
                {
                    posj <- positions[j]
                    typej <- types[j]
                    posdiff <- posj - posi
                    posdiffmod3 <- posdiff %% 3
                    orflength <- posj - posi + 3 # Add in the length of the stop codon
                    if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
                    {
                        # Check if we have already used the stop codon at posj+2 in an ORF
                        numorfs <- length(orfstops)
                        usedstop <- -1
                        if (numorfs > 0)
                        {
                            for (k in 1:numorfs)
                            {
                                orfstopk <- orfstops[k]
                                if (orfstopk == (posj + 2)) { usedstop <- 1 }
                            }
                        }
                        if (usedstop == -1)
                        {
                            orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                            orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                            orflengths <- append(orflengths, orflength, after=length(orflengths))
                        }
                        found <- 1
                        break
                    }
                    if (j == numpositions) { found <- 1 }
                }
            }
        }
    }
    # Sort the final ORFs by start position:
    indices <- order(orfstarts)
    orfstarts <- orfstarts[indices]
    orfstops <- orfstops[indices]
    # Make a plot showing the positions of ORFs in the input sequence:
    # Draw a line at y=0 from 1 to the length of the sequence:
    x <- c(1,nchar(sequence))
    y <- c(0,0)
    plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide # in sequence", ylab="Reading frame", main="Predicted ORFs")
    segments(1,1,nchar(sequence),1)
    segments(1,2,nchar(sequence),2)
    # Add the x-axis at y=0:
    axis(1, pos=0)
    # Add the y-axis labels:
    text(0.9,0.5,"+1")
    text(0.9,1.5,"+2")
    text(0.9,2.5,"+3")
    # Make a plot of the ORFs in the sequence:
    numorfs <- length(orfstarts)
    for (i in 1:numorfs)
    {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        remainder <- (orfstart-1) %% 3
        if    (remainder == 0) # +1 reading frame
        {
            rect(orfstart,0,orfstop,1,col="mediumseagreen",border="navyblue")
        }
        else if (remainder == 1)
        {
            rect(orfstart,1,orfstop,2,col="limegreen",border="navyblue")
        }
        else if (remainder == 2)
        {
            rect(orfstart,2,orfstop,3,col="olivedrab2",border="navyblue")
        }
    }
}

plotORFsinSeqChangeX <- function(sequence, xscale){
    # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
    mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
    orfstarts <- numeric()
    orfstops <- numeric()
    # Make a vector "orflengths" to store the lengths of the ORFs
    orflengths <- numeric()
    # Print out the positions of ORFs in the sequence:
    numpositions <- length(positions) # Find the length of vector "positions"
    # There must be at least one start codon and one stop codon to have an ORF.
    if (numpositions >= 2)
    {
        for (i in 1:(numpositions-1))
        {
            posi <- positions[i]
            typei <- types[i]
            found <- 0
            while (found == 0)
            {
                for (j in (i+1):numpositions)
                {
                    posj <- positions[j]
                    typej <- types[j]
                    posdiff <- posj - posi
                    posdiffmod3 <- posdiff %% 3
                    orflength <- posj - posi + 3 # Add in the length of the stop codon
                    if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
                    {
                        # Check if we have already used the stop codon at posj+2 in an ORF
                        numorfs <- length(orfstops)
                        usedstop <- -1
                        if (numorfs > 0)
                        {
                            for (k in 1:numorfs)
                            {
                                orfstopk <- orfstops[k]
                                if (orfstopk == (posj + 2)) { usedstop <- 1 }
                            }
                        }
                        if (usedstop == -1)
                        {
                            orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                            orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                            orflengths <- append(orflengths, orflength, after=length(orflengths))
                        }
                        found <- 1
                        break
                    }
                    if (j == numpositions) { found <- 1 }
                }
            }
        }
    }
    # Sort the final ORFs by start position:
    indices <- order(orfstarts)
    orfstarts <- orfstarts[indices]
    orfstops <- orfstops[indices]
    # Make a plot showing the positions of ORFs in the input sequence:
    # Draw a line at y=0 from 1 to the length of the sequence:
    x <- c(1,xscale)            # changed it so xscale can be changed manually by input
    y <- c(0,0)
    plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide # in sequence", ylab="Reading frame", main="Predicted ORFs")
    segments(1,1,nchar(sequence),1)
    segments(1,2,nchar(sequence),2)
    # Add the x-axis at y=0:
    axis(1, pos=0)
    # Add the y-axis labels:
    text(0.9,0.5,"+1")
    text(0.9,1.5,"+2")
    text(0.9,2.5,"+3")
    # Make a plot of the ORFs in the sequence:
    numorfs <- length(orfstarts)
    for (i in 1:numorfs)
    {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        remainder <- (orfstart-1) %% 3
        if    (remainder == 0) # +1 reading frame
        {
            rect(orfstart,0,orfstop,1,col="mediumseagreen",border="navyblue")
        }
        else if (remainder == 1)
        {
            rect(orfstart,1,orfstop,2,col="limegreen",border="navyblue")
        }
        else if (remainder == 2)
        {
            rect(orfstart,2,orfstop,3,col="olivedrab2",border="navyblue")
        }
    }
}

LoadGTFAndRemoveLowlyAbundantTranscripts <- function(file , type = "transcript", minTPM = 0, minFPKM = 0, minCov = 0){
    
    # input a loaded stringtie output with transcript annotation
    # input with : read.table(... , sep = "\t", stringsAsFActors = FALSE)
    gtf_input <- read.table(file = file, sep = "\t", stringsAsFactors = FALSE)
    
    
    # type = definition of wether this function shall look at transcripts or exons | transcripts by default
    # TPM, FPKM & Cov  --> cutoffs, per default @ 0
    
    if(type == "transcript"){
        
        tr_keep <- c()
        for(i in c(1:length(rownames(gtf_input)))){
            if(gtf_input$V3[i]=="transcript"){ # for every tr
                if((as.numeric(substr((uncollapse(gtf_input$V9[i])[[1]][6]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][6])-1)))<=minCov){ # check if its cov
                }else if((as.numeric(substr((uncollapse(gtf_input$V9[i])[[1]][8]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][8])-1)))<=minFPKM){ # or FPKM
                }else if((as.numeric(substr((uncollapse(gtf_input$V9[i])[[1]][10]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][10])-1)))<=minTPM){ # or TPM is lower, else add it to tr_keep
                }else{tr_keep <- c(tr_keep,(substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1)))}
            }
        }
        
        rows_keep <- c()
        for(i in c(1:length(rownames(gtf_input)))){
            if(gtf_input$V3[i]=="transcript"){
                if((substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1))%in%tr_keep){rows_keep <- c(rows_keep, i)}
            }
            if(gtf_input$V3[i]=="exon"){
                if((substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1))%in%tr_keep){rows_keep <- c(rows_keep, i)}
            }
            
        }
        
    }
    
    if(type == "exon"){
        if(minCov == 0){
            print("When choosing option Exon, please specify minCov!")
        }
        tr_keep <- c()
        for(i in c(1:length(rownames(gtf_input)))){
            if(gtf_input$V3[i]=="exon"){ # for every tr
                if((as.numeric(substr((uncollapse(gtf_input$V9[i])[[1]][8]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][8])-1)))<=minCov){ # check if its cov
                }else{tr_keep <- c(tr_keep,(substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1)))}
            }
        }
        tr_keep <- unique(tr_keep)
        rows_keep <- c()
        for(i in c(1:length(rownames(gtf_input)))){
            if(gtf_input$V3[i]=="transcript"){
                if((substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1))%in%tr_keep){rows_keep <- c(rows_keep, i)}
            }
            if(gtf_input$V3[i]=="exon"){
                if((substr((uncollapse(gtf_input$V9[i])[[1]][4]), 1, nchar(uncollapse(gtf_input$V9[i])[[1]][4])-1))%in%tr_keep){rows_keep <- c(rows_keep, i)}
            }
            
        }
        
    }
    
    gtf_output <- gtf_input[c(rows_keep),]
    
    
    
    
    
}
# same function as above but gtf does not need to be loaded before

MergeDuplicatesInGtf <- function(gtf_data){  
    
    print("OLD FUNCTION - USE MergeDuplicatesInGtfFreq INSTEAD ")
    # remove every "transcript" row, keep only exon information                    
    exons <- NULL            # make sure 'exons' is empty
    for (i in c(1:(length(rownames(gtf_data))))){   # for every row in gtf_data
        if (grepl("exon",gtf_data$V3[i])==TRUE){  # if the value in the third coloumn == "exon"
            exons <- c(exons, i)                    # add the rowname to the 'exons' string
        }}   
    exons <- gtf_data[c(exons),]                     # 'exons' become a new dataframe including only exon information from gtf_data            
    
    
    # extract all start points
    starts <- unique(exons$V4)
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- exons[exons$V4==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$V5))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the consant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0)))  # dataframe with 3 coloumns: start/end/coverage
    # add rownames later
    for (i in c(1:(length(constant)))){                               # for every constant exon start
        tmp_exon <- exons[exons$V4 == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$V5[x]        # add second location of the exon (should theoretically only be done once, but hey ;) )
            split_v9 <- unlist(strsplit(tmp_exon$V9[x], " "))                # split the last coloumn
            split_v9[8] <- substr(split_v9[8],1,nchar(split_v9[8])-3)       # remove the last character (";")
            const_matrix$V3[i] <- const_matrix$V3[i] + as.numeric(split_v9[8])         # add the coverage value to the start matrix
        }}
    
    
    # make dataframe of exons with variants analog to dataframe of constant exons    
    variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
    
    start_matrix_var <- exons[exons$V4 %in% c(variants),] # extract them
    
    variant_matrix <- NULL
    tmp_var <- NULL
    tmp_unique <- NULL
    for (i in c(1:(length(variants)))){   # for every exon with different ends
        tmp_var <- start_matrix_var[start_matrix_var$V4 == variants[i],]    # extract this dataframe
        tmp_unique <- unique(tmp_var$V5)    # get the unique ends for this start
        for (x in c(tmp_unique)){           # for every differend end
            variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0))) # rbind the start and the differing end to variant_matrix
        }
    }
    
    if ((length(unique(variant_matrix$V2)))==(length(variant_matrix$V2))){       # check if no end is a duplicate 
        print ("MERGING VARIANTS..")
    } else {print ("CAREFUL, SOME DIFFERENT STARTS HAVE THE SAME END")}
    
    
    for (i in c(1:(length(variant_matrix$V2)))){            # for every different end in variant_matrix   (we clarified above in the if statement, that no end is a duplicate of differnet starts)
        tmp_exon <- start_matrix_var[start_matrix_var$V5 == variant_matrix$V2[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){       # for every row
            split_v9 <- unlist(strsplit(tmp_exon$V9[x], " "))                # split the last coloumn
            split_v9[8] <- substr(split_v9[8],1,nchar(split_v9[8])-3)       # remove the last character (";")
            #print(as.numeric(split_v9[8]))                                            # get the coverage value
            variant_matrix$V3[i] <- variant_matrix$V3[i] + as.numeric(split_v9[8])# add the coverage value to the start matrix
        }}
    
    ###  merge to a coverage matrix
    
    colnames(const_matrix) <- c("start","end","cov")
    colnames(variant_matrix) <- c("start","end","cov")
    
    merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
    
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
    
}

MergeDuplicatesInGtfWithFreq <- function(gtf_data){   
    
    # remove every "transcript" row, keep only exon information                    
    exons <- NULL            # make sure 'exons' is empty
    for (i in c(1:(length(rownames(gtf_data))))){   # for every row in gtf_data
        if (grepl("exon",gtf_data$V3[i])==TRUE){  # if the value in the third coloumn == "exon"
            exons <- c(exons, i)                    # add the rowname to the 'exons' string
        }}   
    exons <- gtf_data[c(exons),]                     # 'exons' become a new dataframe including only exon information from gtf_data            
    
    
    # extract all start points
    starts <- unique(exons$V4)
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- exons[exons$V4==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$V5))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the consant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0),c(0)))  # dataframe with 3 coloumns: start/end/coverage/nr.of.transcripts
    # add rownames later
    for (i in c(1:(length(constant)))){                               # for every constant exon start
        tmp_exon <- exons[exons$V4 == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$V5[x]        # add second location of the exon (should theoretically only be done once, but hey ;) )
            split_v9 <- unlist(strsplit(tmp_exon$V9[x], " "))                # split the last coloumn
            split_v9[8] <- substr(split_v9[8],1,nchar(split_v9[8])-3)       # remove the last character (";")
            const_matrix$V3[i] <- const_matrix$V3[i] + as.numeric(split_v9[8])         # add the coverage value to the start matrix
            const_matrix$V4[i] <- const_matrix$V4[i] + 1    # add the number of transcripts
        }}
    
    # check if there are variants
    if(length(starts)==length(constant)){
        ###  merge to a coverage matrix
        
        colnames(const_matrix) <- c("start","end","cov","freq")
        
        merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),])
        
    }else{
        # make dataframe of exons with variants analog to dataframe of constant exons    
        variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
        
        start_matrix_var <- exons[exons$V4 %in% c(variants),] # extract them
        
        variant_matrix <- NULL
        tmp_var <- NULL
        tmp_unique <- NULL
        for (i in c(1:(length(variants)))){   # for every exon with different ends
            tmp_var <- start_matrix_var[start_matrix_var$V4 == variants[i],]    # extract this dataframe
            tmp_unique <- unique(tmp_var$V5)    # get the unique ends for this start
            for (x in c(tmp_unique)){           # for every differend end
                variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0,0))) # rbind the start and the differing end to variant_matrix
            }
        }
        
        if ((length(unique(variant_matrix$V2)))==(length(variant_matrix$V2))){       # check if no end is a duplicate 
            print ("MERGING VARIANTS..")            # if this is a case 
            
            for (i in c(1:(length(variant_matrix$V2)))){            # for every different end in variant_matrix   (we clarified above in the if statement, that no end is a duplicate of differnet starts)
                tmp_exon <- start_matrix_var[start_matrix_var$V5 == variant_matrix$V2[i],]               # extract all rows with this exon 
                for (x in c(1:(length(rownames(tmp_exon))))){       # for every row
                    split_v9 <- unlist(strsplit(tmp_exon$V9[x], " "))                # split the last coloumn
                    split_v9[8] <- substr(split_v9[8],1,nchar(split_v9[8])-3)       # remove the last character (";")
                    #print(as.numeric(split_v9[8]))                                            # get the coverage value
                    variant_matrix$V3[i] <- variant_matrix$V3[i] + as.numeric(split_v9[8])# add the coverage value to the start matrix
                    variant_matrix$V4[i] <- variant_matrix$V4[i] + 1 # add the number of transcripts
                }}
            
        } else {        # if we have different starts which have the same end 
            print ("MERGING VARIANTS ... SOME DIFFERENT STARTS HAVE THE SAME END (FYI)")
            ## extract all starts which have just one end
            #constant_var <- NULL            # make sure 'constant_var' is empty
            #variants_var <- NULL            # make sure 'variants_var' is empty
            #ends <- unique(variant_matrix$V2)# get all unique ends 
            
            #for (i in c(1:(length(ends)))) {  # for every unique exon end in the variant matrix 
            #exon_i_var <- variant_matrix[variant_matrix$V2==ends[i],]   # take all exons with the same end point
            #if (length(unique(exon_i_var$V1))==1){  # if there is only one unique start point
            #     constant_var <- c(constant_var, ends[i])   # add the exon end to the constant_var string
            #  }
            #   else {variants_var <- c(variants_var,ends[i])}               # else: add it to the variants_var string
            #}
            
            
            for (n in c(1:length(rownames(variant_matrix)))){       # for every variant exon
                # print(paste(variant_matrix$V1[i], variant_matrix$V2[i]))    # get exon start and end
                occ_var_i <- start_matrix_var[start_matrix_var$V4==variant_matrix$V1[n] & start_matrix_var$V5 == variant_matrix$V2[n],] # extract all occurences of this exons
                for (y in c(1:(length(rownames(occ_var_i))))){       # for every row
                    split_v9 <- unlist(strsplit(occ_var_i$V9[y], " "))                # split the last coloumn
                    split_v9[8] <- substr(split_v9[8],1,nchar(split_v9[8])-3)       # remove the last character (";")
                    #print(as.numeric(split_v9[8]))                                            # get the coverage value
                    variant_matrix$V3[n] <- variant_matrix$V3[n] + as.numeric(split_v9[8])# add the coverage value to the start matrix
                    variant_matrix$V4[n] <- variant_matrix$V4[n] + 1 # add the number of transcripts
                }}                   }
        
        
        ###  merge to a coverage matrix
        
        colnames(const_matrix) <- c("start","end","cov","freq")
        colnames(variant_matrix) <- c("start","end","cov","freq")
        
        merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
        
    } # this is the closing of the else statement from above when checked if there are variants or not
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
} # now supports the case when there are no variants

MergeGtfInformation <- function(merge_list, replications){
    
    print("OLD FUNCTION - USE MergeGTFInformationWithFreq INSTEAD ")
    
    # This function bases on the Other MergeGtf - function by me
    # merge_list has to be made in the following way:
    # merge_list <- rbind(chow_1_merged.gtf[1:35,], chow_2_merged.gtf[1:33,], chow_3_merged.gtf[1:35,])
    # extract all start points
    starts <- unique(merge_list$start)
    starts
    
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- merge_list[merge_list$start==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$end))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the consant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0)))  # dataframe with 3 coloumns: start/end/coverage
    # add rownames later
    for (i in c(1:(length(constant)))){                               # for every constant exon start
        tmp_exon <- merge_list[merge_list$start == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$end[x]                             # add second location of the exon (should theoretically only be done once, but hey ;) )
            const_matrix$V3[i] <- const_matrix$V3[i] + tmp_exon$cov[x]        # add the coverage value to the start matrix
        }}
    
    
    # make dataframe of exons with variants analog to dataframe of constant exons    
    variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
    
    start_matrix_var <- merge_list[merge_list$start %in% c(variants),] # extract them
    
    variant_matrix <- NULL
    tmp_var <- NULL
    tmp_unique <- NULL
    for (i in c(1:(length(variants)))){   # for every exon with different ends
        tmp_var <- start_matrix_var[start_matrix_var$start == variants[i],]    # extract this dataframe
        tmp_unique <- unique(tmp_var$end)    # get the unique ends for this start
        for (x in c(tmp_unique)){           # for every differend end
            variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0))) # rbind the start and the differing end to variant_matrix
        }
    }
    
    if ((length(unique(variant_matrix$V2)))==(length(variant_matrix$V2))){       # check if no end is a duplicate 
        print ("MERGING VARIANTS..")
    } else {print ("CAREFUL, SOME DIFFERENT STARTS HAVE THE SAME END")}
    
    
    for (i in c(1:(length(variant_matrix$V2)))){            # for every different end in variant_matrix   (we clarified above in the if statement, that no end is a duplicate of differnet starts)
        tmp_exon <- start_matrix_var[start_matrix_var$end == variant_matrix$V2[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){       # for every row
            variant_matrix$V3[i] <- variant_matrix$V3[i] + tmp_exon[x,3] # add the coverage value to the start matrix
        }}
    
    ###  merge to a coverage matrix
    
    colnames(const_matrix) <- c("start","end","cov")
    colnames(variant_matrix) <- c("start","end","cov")
    
    merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
    
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # get the mean coverage
    
    sorted_merged.gtf$cov <- sorted_merged.gtf$cov/as.integer(replications)
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
    
}

MergeGtfInformationWithFreq <- function(merge_list, replications){
    
    # This function bases on the Other MergeGtf - function by me
    # merge_list has to be made in the following way:
    # merge_list <- rbind(chow_1_merged.gtf[1:35,], chow_2_merged.gtf[1:33,], chow_3_merged.gtf[1:35,])
    # extract all start points
    starts <- unique(merge_list$start)
    starts
    
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- merge_list[merge_list$start==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$end))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the consant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0),c(0)))  # dataframe with 3 coloumns: start/end/coverage/tr.frequency
    # add rownames later
    
    for (i in c(1:(length(constant)))){                               # for every constant exon start
        tmp_exon <- merge_list[merge_list$start == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$end[x]                             # add second location of the exon (should theoretically only be done once, but hey ;) )
            const_matrix$V3[i] <- const_matrix$V3[i] + tmp_exon$cov[x]          # add the coverage value to the start matrix
            const_matrix$V4[i] <- const_matrix$V4[i] + tmp_exon$freq[x]     # add frequency
        }}
    
    
    
    
    # make dataframe of exons with variants analog to dataframe of constant exons    
    variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
    
    
    
    start_matrix_var <- merge_list[merge_list$start %in% c(variants),] # extract them
    
    variant_matrix <- NULL
    tmp_var <- NULL
    tmp_unique <- NULL
    
    #
    
    for (i in c(1:(length(variants)))){   # for every exon with different ends
        tmp_var <- start_matrix_var[start_matrix_var$start == variants[i],]    # extract this dataframe
        tmp_unique <- unique(tmp_var$end)    # get the unique ends for this start
        for (x in c(tmp_unique)){           # for every differend end
            variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0,0))) # rbind the start and the differing end to variant_matrix
        }
    }
    
    
    for (i in c(1:(length(rownames(variant_matrix))))){ # for every variant
        # print(paste(variant_matrix$V1[i], variant_matrix$V2[i]))    # get exon start and end
        occ_var_i <- start_matrix_var[start_matrix_var$start==variant_matrix$V1[i] & start_matrix_var$end == variant_matrix$V2[i],] # extract all rows with this exon
        for (n in c(1:(length(rownames(occ_var_i))))){      # for every row
            variant_matrix$V3[i] <- variant_matrix$V3[i] + occ_var_i[n,3] # add coverage to variant matrix
            variant_matrix$V4[i] <- variant_matrix$V4[i] + occ_var_i[n,4] # add freq
            
        }
    }
    
    
    
    ###  merge to a coverage matrix
    
    colnames(const_matrix) <- c("start","end","cov","freq")
    colnames(variant_matrix) <- c("start","end","cov","freq")
    
    merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
    
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # get the mean coverage
    
    sorted_merged.gtf$cov <- sorted_merged.gtf$cov/as.integer(replications)
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
    
}

MergeGtfInformationWithFreqInclCondition <- function(merge_list, replications){
    
    # replications has to have at least 3 values c(rep_all, rep_cond1, rep_cond2)
    # this is an extention of MergeGTFInformationWithFreq
    # merge_list has to have an additional col specifying the condition (right now works with 2 conditions)
    
    conditions <- as.character(unique(merge_list$cond))
    
    # extract all start points
    starts <- unique(merge_list$start)
    starts
    
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- merge_list[merge_list$start==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$end))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the constant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0), c(0), c(0),c(0)))  # dataframe with 6 coloumns: start/end/coverage//tr.frequency/cond1_cov/cond2_cov
    
    
    # add rownames later
    for (i in c(1:(length(constant)))){ # for every constant exon start
        tmp_exon <- merge_list[merge_list$start == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$end[x]                             # add second location of the exon (should theoretically only be done once, but hey ;) )
            const_matrix$V3[i] <- const_matrix$V3[i] + tmp_exon$cov[x]          # add the coverage value to the start matrix
            const_matrix$V4[i] <- const_matrix$V4[i] + tmp_exon$freq[x]     # add frequency
            if(tmp_exon$cond[x]==conditions[1]){const_matrix$V5[i] <- const_matrix$V5[i] + tmp_exon$cov[x]} # add freq if it is cond 1
            if(tmp_exon$cond[x]==conditions[2]){const_matrix$V6[i] <- const_matrix$V6[i] + tmp_exon$cov[x]} # add freq if it is cond 2
            
        }}
    
    
    
    
    # make dataframe of exons with variants analog to dataframe of constant exons    
    
    variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
    
    
    
    start_matrix_var <- merge_list[merge_list$start %in% c(variants),] # extract them
    
    variant_matrix <- NULL
    tmp_var <- NULL
    tmp_unique <- NULL
    
    #
    
    for (i in c(1:(length(variants)))){   # for every exon with different ends
        tmp_var <- start_matrix_var[start_matrix_var$start == variants[i],]    # extract this dataframe
        tmp_unique <- unique(tmp_var$end)    # get the unique ends for this start
        for (x in c(tmp_unique)){           # for every differend end
            variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0,0,0,0))) # rbind the start and the differing end to variant_matrix
        }
    }
    
    
    for (i in c(1:(length(rownames(variant_matrix))))){ # for every variant
        # print(paste(variant_matrix$V1[i], variant_matrix$V2[i]))    # get exon start and end
        
        occ_var_i <- start_matrix_var[start_matrix_var$start==variant_matrix$V1[i] & start_matrix_var$end == variant_matrix$V2[i],] # extract all rows with this exon
        for (n in c(1:(length(rownames(occ_var_i))))){      # for every row
            
            variant_matrix$V3[i] <- variant_matrix$V3[i] + occ_var_i[n,3] # add coverage to variant matrix
            variant_matrix$V4[i] <- variant_matrix$V4[i] + occ_var_i[n,4] # add freq
            if(occ_var_i$cond[n]==conditions[1]){variant_matrix$V5[i] <- variant_matrix$V5[i] + occ_var_i$cov[n]} # add freq if it is cond 1
            if(occ_var_i$cond[n]==conditions[2]){variant_matrix$V6[i] <- variant_matrix$V6[i] + occ_var_i$cov[n]} # add freq if it is cond 2
            
        }
    }
    
    
    
    ###  merge to a coverage matrix
    
    colnames(const_matrix) <- c("start","end","cov","freq", "cov_cond1", "cov_cond2")
    colnames(variant_matrix) <- c("start","end","cov","freq", "cov_cond1", "cov_cond2")
    
    merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
    
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # get the mean coverage
    sorted_merged.gtf$cov <- sorted_merged.gtf$cov/as.integer(replications[1])
    sorted_merged.gtf$cov_cond1 <- sorted_merged.gtf$cov_cond1/as.integer(replications[2])
    sorted_merged.gtf$cov_cond2 <- sorted_merged.gtf$cov_cond2/as.integer(replications[3])
    
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
    
}

ExportExonDataFrameAsGtf <- function(exon_df, chr, strandedness, gene_id, file){
    #exon_df = dataframe as obtained by my other function MergeDuplicatesInGtf or MergeGtfInformation
    # chr = chromosome
    # strandedness = "+" / "-"
    # gene_id = Name of the Gene
    # file = desired name of the output file, dont forget '.gtf'
    
    #make dataframe with 9 coloumns
    
    row_len <- c(1:length(rownames(exon_df)))   # same number of rows as the original dataframe
    export.gtf <- as.data.frame(cbind(row_len,row_len,row_len,row_len,row_len,row_len,row_len,row_len,row_len))
    
    export.gtf[,1] <- chr
    export.gtf[,2] <- c("SmbAK")   # Stringtie merged by Alex Knierim ;)
    export.gtf[,3] <- c("exon")
    export.gtf[,4] <- c(exon_df$start)
    export.gtf[,5] <- c(exon_df$end)
    export.gtf[,6] <- c(1000) # don't know what that means
    export.gtf[,7] <- c(strandedness) # strandedness
    export.gtf[,8] <- c(".")
    
    # check if input has frequency information
    
    if (length(colnames(exon_df))>=4) {     # if there are 4 coloumns (last coloumn is frequency information)
        print("Including Frequency Information")
        for (i in row_len){
            c9 <- paste('gene_id "',gene_id, '";transcript_id ',paste('"TR',i,'"',sep = ""),';exon_number ',as.character(i),'; ','cov "', exon_df$cov[i],'"; ','freq "', exon_df$freq[i],'"; ', sep = "") # 
            export.gtf[i,9] <- as.character(c9)
        }                    
    }
    else{               # if there are 3 coloumns (no frequency informaiton)
        print("No Frequency Information Given, Proceed Without ..")
        # ninth coloumn is tricky
        for (i in row_len){
            c9 <- paste('gene_id "',gene_id, '";transcript_id ',paste('"TR',i,'"',sep = ""),';exon_number ',as.character(i),'; ','cov "', exon_df$cov[i],'"; ', sep = "") # 
            export.gtf[i,9] <- as.character(c9)
        }                    
    }
    
    #export
    write.table(export.gtf, file = file, sep = "\t", col.names = F, row.names = F, quote = F)
}

ExportExonDataFrameAsGtfInclCondition <- function(exon_df, chr, strandedness, gene_id, file){
    
    #exon_df = dataframe as obtained by my other function MergeDuplicatesInGtfIncludingCondition
    # chr = chromosome
    # strandedness = "+" / "-"
    # gene_id = Name of the Gene
    # file = desired name of the output file, dont forget '.gtf'
    
    #make dataframe with 9 coloumns
    
    row_len <- c(1:length(rownames(exon_df)))   # same number of rows as the original dataframe
    export.gtf <- as.data.frame(cbind(row_len,row_len,row_len,row_len,row_len,row_len,row_len,row_len,row_len))
    
    export.gtf[,1] <- chr
    export.gtf[,2] <- c("SmbAK")   # Stringtie merged by Alex Knierim ;)
    export.gtf[,3] <- c("exon")
    export.gtf[,4] <- c(exon_df$start)
    export.gtf[,5] <- c(exon_df$end)
    export.gtf[,6] <- c(1000) # don't know what that means
    export.gtf[,7] <- c(strandedness) # strandedness
    export.gtf[,8] <- c(".")
    
    # 9th coloumn
    
    for (i in row_len){
        c9 <- paste('gene_id "',gene_id, '";transcript_id ',paste('"TR',i,'"',sep = ""),';exon_number ',as.character(i),'; ','cov "', exon_df$cov[i],'"; ','freq "', exon_df$freq[i],'"; ','cov_cond1 "',exon_df$cov_cond1[i],'"; ','cov_cond1 "',exon_df$cov_cond2[i],'";', sep = "") # 
        export.gtf[i,9] <- as.character(c9)
    }
    
    #export
    write.table(export.gtf, file = file, sep = "\t", col.names = F, row.names = F, quote = F)
}

MergeGtfInformationWithFreqInclCondition <- function(merge_list, replications){
    
    # replications has to have at least 3 values c(rep_all, rep_cond1, rep_cond2)
    # this is an extention of MergeGTFInformationWithFreq
    # merge_list has to have an additional col specifying the condition (right now works with 2 conditions)
    
    conditions <- as.character(unique(merge_list$cond))
    
    # extract all start points
    starts <- unique(merge_list$start)
    starts
    
    # find out wether an exon is constant (allways the same start/end location) or has variants (differing end points)
    # This is just a computational approach: In reality an exon can have different begins and the same end. Since the start location
    # is used for exon definition, different start positions will make different exons and not variants of the same exons. 
    constant <- NULL            # make sure 'constant' is empty
    variants <- NULL            # make sure 'variants' is empty
    
    for (i in c(1:(length(starts)))) {  # for every unique exon start
        exon_i <- merge_list[merge_list$start==starts[i],]   # take all exons with the same start point
        if (length(unique(exon_i$end))==1){  # if there is only one unique end point
            constant <- c(constant,i)   # add the number of the exon to the constant string
        }
        else {variants <- c(variants,i)}               # else: add it to the variants string
    }
    
    # make dataframe of constant exons with start/end & combined coverage
    constant <- starts[constant]                # overwrite the number of the exons with the constants (start locations)
    const_matrix <- as.data.frame(cbind(constant,c(0),c(0), c(0), c(0),c(0)))  # dataframe with 6 coloumns: start/end/coverage//tr.frequency/cond1_cov/cond2_cov
    
    
    # add rownames later
    for (i in c(1:(length(constant)))){ # for every constant exon start
        tmp_exon <- merge_list[merge_list$start == constant[i],]               # extract all rows with this exon 
        for (x in c(1:(length(rownames(tmp_exon))))){   # for every row
            const_matrix$V2[i] <- tmp_exon$end[x]                             # add second location of the exon (should theoretically only be done once, but hey ;) )
            const_matrix$V3[i] <- const_matrix$V3[i] + tmp_exon$cov[x]          # add the coverage value to the start matrix
            const_matrix$V4[i] <- const_matrix$V4[i] + tmp_exon$freq[x]     # add frequency
            if(tmp_exon$cond[x]==conditions[1]){const_matrix$V5[i] <- const_matrix$V5[i] + tmp_exon$cov[x]} # add freq if it is cond 1
            if(tmp_exon$cond[x]==conditions[2]){const_matrix$V6[i] <- const_matrix$V6[i] + tmp_exon$cov[x]} # add freq if it is cond 2
            
        }}
    
    
    
    
    # make dataframe of exons with variants analog to dataframe of constant exons    
    
    variants <- starts[variants]        # overwrite the number of the exons with the variants (start locations)
    
    
    
    start_matrix_var <- merge_list[merge_list$start %in% c(variants),] # extract them
    
    variant_matrix <- NULL
    tmp_var <- NULL
    tmp_unique <- NULL
    
    #
    
    for (i in c(1:(length(variants)))){   # for every exon with different ends
        tmp_var <- start_matrix_var[start_matrix_var$start == variants[i],]    # extract this dataframe
        tmp_unique <- unique(tmp_var$end)    # get the unique ends for this start
        for (x in c(tmp_unique)){           # for every differend end
            variant_matrix <- as.data.frame(rbind(variant_matrix, c(variants[i],x,0,0,0,0))) # rbind the start and the differing end to variant_matrix
        }
    }
    
    
    for (i in c(1:(length(rownames(variant_matrix))))){ # for every variant
        # print(paste(variant_matrix$V1[i], variant_matrix$V2[i]))    # get exon start and end
        
        occ_var_i <- start_matrix_var[start_matrix_var$start==variant_matrix$V1[i] & start_matrix_var$end == variant_matrix$V2[i],] # extract all rows with this exon
        for (n in c(1:(length(rownames(occ_var_i))))){      # for every row
            
            variant_matrix$V3[i] <- variant_matrix$V3[i] + occ_var_i[n,3] # add coverage to variant matrix
            variant_matrix$V4[i] <- variant_matrix$V4[i] + occ_var_i[n,4] # add freq
            if(occ_var_i$cond[n]==conditions[1]){variant_matrix$V5[i] <- variant_matrix$V5[i] + occ_var_i$cov[n]} # add freq if it is cond 1
            if(occ_var_i$cond[n]==conditions[2]){variant_matrix$V6[i] <- variant_matrix$V6[i] + occ_var_i$cov[n]} # add freq if it is cond 2
            
        }
    }
    
    
    
    ###  merge to a coverage matrix
    
    colnames(const_matrix) <- c("start","end","cov","freq", "cov_cond1", "cov_cond2")
    colnames(variant_matrix) <- c("start","end","cov","freq", "cov_cond1", "cov_cond2")
    
    merged.gtf <- rbind(const_matrix[1:c(length(rownames(const_matrix))),], variant_matrix[1:c(length(rownames(variant_matrix))),])
    
    # sort the coverage matrix
    
    sorted_merged.gtf <- merged.gtf[order(merged.gtf$start, merged.gtf$end),]
    
    # get the mean coverage
    sorted_merged.gtf$cov <- sorted_merged.gtf$cov/as.integer(replications[1])
    sorted_merged.gtf$cov_cond1 <- sorted_merged.gtf$cov_cond1/as.integer(replications[2])
    sorted_merged.gtf$cov_cond2 <- sorted_merged.gtf$cov_cond2/as.integer(replications[3])
    
    
    # output the coverage matrix
    
    return(sorted_merged.gtf)
    
}

AppendCovAndUniqueIdentifierAndExport <- function(gtf,meta_data,filename){
    for (i in c(1:(length(rownames(gtf))))){  # for every rownames
        if (gtf$V3[i] == "exon"){           # for each exon
            for (x in c(1:(length(rownames(meta_data))))){   # loop over the whole identifier matrix
                if (gtf$V4[i] == meta_data$start[x]){        # if start matches
                    if (gtf$V5[i] == meta_data$end[x]){      # if end matches
                        gtf[i,10] <- x                      # add identifier
                        gtf[i,11] <- meta_data$cov[x]
                    } 
                }
            }       
        }
        else gtf[i,10] <- NA    
    } 
    # append the identifier & the coverage to the 9th coloumn
    for (i in c(1:(length(rownames(gtf))))){
        if (gtf$V3[i] == "exon"){
            gtf$V9[i] <- paste(gtf$V9[i],' uniq_id "',as.character(gtf[i,10]),'"; ','mean_cov "',as.character(gtf[i,11]),'";', sep = "")
        }           
    }
    # remove duplicate information
    gtf_out <- gtf[,1:9]
    # export 
    write.table(gtf_out, file = filename, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

ConvertGTFtoGFFandExport <- function(transcripts.gtf, meta_data, filename){
    
    # meta_data = unique exon identifier with cols: start/end/cov/freq [only the first two are being used atm]
    # transcripts.gtf = concatenated transcript files of stringtie output, CAREFUL: There should be no duplicate transcript id
    # filename = desired output
    
    # gff file has following coloumns 
    # chr origin.script aggregate_gene/exonic_part  start   stop    .   +   .   (transcripts "xx.1+xx.2"); gene_id "xx"
    # create gff template
    
    col <- c(1:(length(rownames(meta_data))+1))
    gff <- as.data.frame(cbind(col,col,col,col,col,col,col,col,col))
    colnames(gff) <- c("chr","origin","structure","start","stop","pt1","strand","pt2","tr_info")
    
    # consistent cols
    
    gff$chr <- transcripts.gtf$V1[1]
    gff$origin <- "RFunctionAK"  #RFunction Alex Knierim ;)
    gff$structure <- c("exonic_part")
    gff$structure[1] <- c("aggregate_gene")
    gff$pt1 <- c(".")
    gff$pt2 <- c(".")
    gff$strand <- transcripts.gtf$V7[1]
    
    # first row first (important for rereferencing which gene the transcripts belong to)
    
    gff$start[1] <- meta_data[1,1]
    gff$stop[1] <- meta_data[c(length(meta_data[,2])),2]
    
    gene <- unlist(strsplit(transcript.gtf$V9[1], " "))[2]                # split the last coloumn
    gene <- substr(gene,1,nchar(gene)-1)
    gff$tr_info[1] <- paste('gene_id "',as.character(gene),'"')
    
    # start and stops of exons
    
    gff$start[2:(length(rownames(meta_data))+1)] <- meta_data[,1]
    gff$stop[2:(length(rownames(meta_data))+1)] <- meta_data[,2]
    
    # end of last coloumn for exons exonic_part_number "x"; gene_id "gene"
    
    for (i in c(1:(length(rownames(meta_data))))){
        gff$tr_info[(i+1)] <-paste('exonic_part_number "',as.character(as.numeric(i)),'"; gene_id "',as.character(gene),'"', sep="")
    }
    
    ###
    
    #for every beginn after row 2 
    #   loop over every row of transcripts and if start and end match, add transcript identifier 
    #  add list of transcript identifiers to last coloumn
    
    
    
    for (i in 2:(length(rownames(gff)))){       # for every exon in the gff file
        exon_tr <- transcripts.gtf[transcripts.gtf$V4 == gff$start[i] & transcripts.gtf$V5 == gff$stop[i] & transcripts.gtf$V3 == "exon",] # extract every occurence of this exon
        trs_tmp <- NULL    
        for (x in c(1:length(rownames(exon_tr)))){  # loop over every line and extract TR identifier in a seperate variable
            tr_tmp <- unlist(strsplit(exon_tr$V9[x], " "))[4] #extract the single transcript                                    
            tr_tmp <- substr(tr_tmp,1,nchar(tr_tmp)-1)  # remove the last semicolon
            trs_tmp <- c(trs_tmp, as.character (tr_tmp)) # add to list
        }   
        trans_tmp <- paste('transcripts "') # start transcript list suited for gff format
        for (n in c(1:length(trs_tmp))){   # for every transcript
            if (n == length(trs_tmp)){  # last should be added without "+"
                trans_tmp <- paste(trans_tmp,'',as.character(trs_tmp[n]),'";', sep = "")
            }
            else {trans_tmp <- paste(trans_tmp,'',as.character(trs_tmp[n]),'+', sep = "")  # add it with a plus sign
            }}
        gff$tr_info[i] <- paste(trans_tmp,gff$tr_info[i],sep = "")
    }
    
    # export it
    
    write.table(file = filename, x = gff, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
}

plotMotifsInSequence <- function(sequence, motif_table, return_table = FALSE){
    # motif_table should have following coloumns: Sequence (Upper case), Name, Condition (tab, seperated), there is a maximum of 5
    # conditions supported, and the motifs have to be in the same order as the conditions (first the condition1 motifs, then the condition2 motifs etc.)
    # plot basic graph
    x  <- c(1,nchar(sequence)+200)  # more place for legend later
    y <- c(0,0)
    plot(x, y, ylim=c(0,1), type = "l", axes=FALSE, xlab="Amino Acids # in sequence", ylab = "Sequence motifs",
         main="Included Sequence motifs in Protein Sequence")
    axis(1, pos=0)
    
    # color palettes
    # get colour palettes (for each of which a different color palette is needed) (5 are supported)
    # colfunc <- colorRampPalette(c("darkblue", "deepskyblue"))
    
    col_cond_1 <- colorRampPalette(c("darkblue", "darkturquoise"))( length(rownames(motif_table[motif_table$V3==1,])))
    col_cond_2 <-  colorRampPalette(c("forestgreen","lawngreen"))( length(rownames(motif_table[motif_table$V3==2,])))
    col_cond_3 <- colorRampPalette(c("red1","red4"))( length(rownames(motif_table[motif_table$V3==3,])))
    col_cond_4 <- colorRampPalette(c("yellow4","yellow1"))( length(rownames(motif_table[motif_table$V3==4,])))
    col_cond_5 <- colorRampPalette(c("grey0","grey80"))( length(rownames(motif_table[motif_table$V3==5,])))
    col_conds <- c(col_cond_1, col_cond_2, col_cond_3, col_cond_4, col_cond_5)
    
    
    if(length(unique(motif_table$V3))>=6){
        print("More than 5 conditions are not supported.")
    }
    
    # loop over every motif (also possible if there is only one) & get information & plot it
    included <- c()
    for (i in c(1:length(rownames(motif_table)))){  # for every motif
        motif <- motif_table$V1[i]
        occurrences <- matchPattern(motif, sequence) # find every occurrence(s) in the sequence
        motifstarts <- start(occurrences)
        motifends   <- end(occurrences)
        
        numoccurrences <- length(occurrences) # find total number of occurence of this motiv in the sequence
        
        if (numoccurrences == 0){           # if the sequence is not found
            #print(paste("Sequence ",motif_table$V2[i]," not found... Are both strings upper case/lower case?"))
            included <- c(included, "no")    
        }
        else{
            
            mot_col <- col_conds[i]
            
            rect(motifstarts,0,motifends,1,col=mot_col)
            included <- c(included, "yes")
        }
        
    }
    legend(legend =c(motif_table$V2), x = nchar(sequence)+10, y = 1,#lty=c(2,2),
           lwd=c(5,5),cex=0.6,col=col_conds[1:length(rownames(motif_table))], bty = "n" ) # gives the legend lines the correct color and width
    if(return_table == TRUE){
        motif_table <- cbind(motif_table, included)
        return(motif_table)
    }
    
}
# removed verbose output and shrinked legend

GetORFAndProteinSeqFromConcatGTFAndMotifFinder <- function(exons_fasta, conc_tr.gtf, pdfname, orfseqname, protseqname, motif_table){
    
    ### prepare exon information (sequences)
    
    exons_fasta <- read.table(file = exons_fasta, stringsAsFactors = FALSE) # load fasta file
    exons <- exons_fasta$V1[(seq(1, (length(rownames(exons_fasta))-1), 2))] # get first row (exon names, SHOULD BE NUMERIC AND COHERENT WITH THE UNIQ IDs)
    exons_num <- c()
    for (i in exons){ # remove the > & the factor   
        exons_num <- c(exons_num, as.numeric(paste0(strsplit(i, '')[[1]][-1], collapse = '')))
    }
    fastas <- as.character(exons_fasta$V1[(seq(2,length(rownames(exons_fasta)),2))]) # get fastasequences
    
    exons_fasta <- as.data.frame(rbind(exons_num, fastas) ) # merge first and second row together
    colnames(exons_fasta) <- exons_num
    
    ## since I am unable to remove the factors easily, I export and reimport it with stringsAsFactors (hopefully noone else ever reads this)
    
    #write.table(file = "xxx_tmp_pls_delete_xxx", exons_fasta, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    #exons_fasta <- read.table(file = "xxx_tmp_pls_delete_xxx", sep = "\t", stringsAsFactors = FALSE)
    
    # new variants
    
    exons_fasta <- as.data.frame(taRifx::remove.factors(exons_fasta))
    
    
    ### prepare transcript table
    tr.gtf <- read.table(file = conc_tr.gtf, sep = "\t", stringsAsFactors = FALSE) # load concat tr information
    
    # check if it is a reverse gene
    if(unique(tr.gtf$V7) == "-"){ # if strandedness in input is "-"
        reverse_gene <- as.character("yes") # set rev = TRUE
    }else if(unique(tr.gtf$V7)=="+"){
        reverse_gene <- as.character("no")#everything is fine
    }else{print("CAREFUL; DIFFERENT STRANDEDNESS FOR DIFFERENT TRANSCRIPT")}
    
    
    tr_names <- c() 
    for (i in c(1:length(rownames(tr.gtf)))){      # get transcript names
        if (tr.gtf$V3[i]=="transcript"){
            tr_split <- (strsplit(tr.gtf$V9[i], ";"))[[1]][2]   # get tr name
            tr_names <- c(tr_names, gsub(' transcript_id ','',tr_split))
        }
    }
    
    tr <- data.frame(tr_names,c(NA),c(NA),c(NA),c(NA), stringsAsFactors = FALSE) # make new dataframe
    colnames(tr) <- c("V1","V2","V3","V4","V5")
    
    ### merge transcript sequences toegether
    
    for (i in c(1:length(rownames(tr)))){ # for every transcript
        uniqs <- c()
        for (x in (c(1:length(rownames(tr.gtf))))) {# for every possible row (in tr.gtf)
            if(tr.gtf$V3[x]=="exon"){ # for every exon
                tr_split <- gsub(' transcript_id ','',((strsplit(tr.gtf$V9[x],";"))[[1]][2]))  # get belonging transcript    
                if(tr_split == tr$V1[i]){ #if it is the same --> ad uniq identifier
                    uniq_split <- (strsplit(tr.gtf$V9[x], ";"))[[1]][5]
                    #uniqs <- c(uniqs, as.character(gsub(' uniq_id ','', uniq_split)))
                    uniqs <- as.character(paste(uniqs, as.character(gsub(' uniq_id ','', uniq_split)), sep = " "))
                }   
            }                      
        }
        tr$V2[i] <- uniqs       # add exon sequence as input for second coloumn
    }
    
    # add sequence, orf, translation & plots
    orfless <- c()
    for (i in c(1:(length(rownames(tr))))){     # for every transcript
        
        # in the beginning, open pdf
        if(i ==1){
            pdf(file = pdfname)
            par( mfrow = c(2,1))
        }
        
        exons_i <- as.character((tr[i,2]))  # extract the composition of exons
        exons_i <- (strsplit(exons_i, " "))
        exons_i <- as.numeric(exons_i[[1]])
        seq <- c()
        
        for (x in exons_i){                         # for every exon in the transcript
            #print(as.character(exons_fasta[2,x]))        
            seq_x <- as.character(exons_fasta[2,x]) # get sequence
            seq <- c(seq, seq_x)                # add sequence           
        }
        seq <- paste(seq, collapse="")     # collapse
        
        
        if (reverse_gene == "yes"){ # if it is a minus strand
            seq <- DNAString(seq)
            seq <- reverseComplement(seq) # make reverse complement
        }
        
        tr$V3[i] <- tolower(as.character(seq))
        # if no ORF is found
        
        if (is.na(findORFsinSeq(tr$V3[i])[[3]][1])){
            tr$V4[i] <- as.character("NO ORF FOUND")
            tr$V5[i] <- as.character("NO ORF FOUND")
            orfless <- c(orfless,i)
        }else{
            
            orf_dt <- as.data.frame(findORFsinSeq(tr$V3[i]))    # get possible ORFS
            colnames(orf_dt) <- c("start", "stop", "length")    # rename the coloumns
            orf_dt <- orf_dt[order(orf_dt$length, decreasing = TRUE),]  # order by decreasing
            orf_start <- orf_dt$start[1]      # longest is in the first rows
            orf_end <- orf_dt$stop[1]
            tr$V4[i] <- c2s((s2c(tr$V3[i]))[orf_start:orf_end])
            tr$V5[i]<- c2s(seqinr::translate(s2c(tr$V4[i])))
            
        }
        
        if (i == length(rownames(tr))){   # after the last loop
            
            for (n in c(1:(length(rownames(tr))))){   #plot everything
                
                if (n %in% c(orfless)){
                    # don't plot anything
                }else{
                    x_max <- max(nchar(tr$V3))# get longest sequence
                    #plotPotentialStartsAndStopsChangeX(tr$V3[n], xscale = x_max)
                    #title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                    plotORFsinSeqChangeX(tr$V3[n], xscale = x_max)
                    title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                    plotMotifsInSequence(tr$V5[n], motif_table = motif_table)
                    title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                }
            }
            
            
            dev.off()
            output <- c(tr$V4)           # output orfs   
            output <- as.list(output)
            names(output) <- c(as.character(tr$V1))
            write.fasta(sequences = output, names = names(output), file.out = orfseqname)
            
            output <- c(tr$V5)              # output AA sequence
            output <- as.list(output)
            names(output) <- c(as.character(tr$V1))
            
            write.fasta(sequences = output, names = names(output), file.out = protseqname)
        }
        
    }
    
    # add peptide mass
    tr <- cbind(tr,NA)
    colnames(tr)[6]<-"V6"
    
    for(i in c(1:length(rownames(tr)))){
        tr$V6[i]<-CalculatePeptideMass(peptide_sequence = tr$V5[i], unit = "kDa")
    }
    
    # add transcript abundance
    tr <- cbind(tr,NA)
    colnames(tr)[7]<-"V7"
    
    trans_data <- tr.gtf[tr.gtf$V3 == "transcript",]
    
    
    for(i in c(1:length(rownames(tr)))){
        trans_name <- tr$V1[i]
        for(n in allrows(trans_data)){
            tr_tmp <- substr(uncollapse(trans_data$V9[n])[[1]][4], 1, nchar(uncollapse(trans_data$V9[n])[[1]][4])-1)
            if (tr_tmp == trans_name){
                tr$V7[i] <- paste(uncollapse(trans_data$V9[n])[[1]][5:10],collapse=" ")
            }
        }
        
    }
    
    # add in the 8th coloumn: the orf that starts with the signal peptide
    s_pept <- NA
    
    for(i in allrows(motif_table)){
        if(motif_table$V2[i]%in%c("s-pept","S-Pept","signal peptide","spept","Signalpeptid","s pept","s.-pept.","S.-Pept.","S-Peptide","S-Peptid","S.-Peptid","S.-peptide")){
            s_pept<- motif_table$V1[i]
        }
    }
    
    if(is.na(s_pept)){
        print("no signal peptide in motif input")
        # add nothing to 8th coloumkn
        tr <- taRifx::remove.factors(as.data.frame(cbind(tr,"ANALYSIS NOT POSSIBLE")))
        
    }else{ # look @ every ORF --> does it have the S-Pept? --> if yes --> put it into 8th coloumn, if not look for other orfs
        
        col_8 <- c()    
        
        for(i in allrows(tr)){
            #print(i)
            #i <- 1
            if(grepl(pattern = s_pept, x = tr$V5[i])){ # if s-pept is already in there
                col_8 <- c(col_8, tr$V5[i])
            }else{
                orf_dt <- as.data.frame(findORFsinSeq(tr$V3[i]))    # get possible ORFS
                colnames(orf_dt) <- c("start", "stop", "length")    # rename the coloumns
                orf_dt <- orf_dt[order(orf_dt$length, decreasing = TRUE),]  # order by decreasing
                orf_dt <-orf_dt[orf_dt$length>=TripleLength(s_pept),] # remove every orf shorter than 3x s-pept
                orf_dt <- orf_dt[-1,] # remove longest (we already used it)
                if(nrow(orf_dt)<1){
                    col_8 <- c(col_8,"-")
                }else{
                    s_pept_found <- FALSE
                    n <- 0
                    while(s_pept_found==FALSE){
                        n <- n + 1
                        orf_start <- orf_dt$start[n]      
                        orf_end <- orf_dt$stop[n]
                        coding_seq <- c2s((s2c(tr$V3[i]))[orf_start:orf_end])
                        pept_seq <- c2s(seqinr::translate(s2c(coding_seq)))
                        if(grepl(pattern = s_pept, x = pept_seq)){
                            col_8 <- c(col_8, pept_seq)
                            tr$V1[i] <- paste(tr$V1[i],"_(mo)",collapse = "", sep  ="")
                            tr$V4[i] <- paste(tr$V4[i],coding_seq, paste=" ", collapse=" ")
                            s_pept_found <- "yes"
                            print("yes")
                        }else if(n == nrow(orf_dt)){ # if the last row is reached and still no orf with the s-pept is found -->
                            col_8 <- c(col_8,"-")
                            s_pept_found <- "no"
                            print("no")
                        }
                    }
                    
                }
            }
            
        }
        
        tr <- taRifx::remove.factors(as.data.frame(cbind(tr,col_8)))
        
        
    }
    
    
    colnames(tr)[8] <- "V8"
    return(tr)
}
# added possibility of two orfs

GetORFAndProteinSeqFromConcatGTFAndMotifFinder <- function(exons_fasta, conc_tr.gtf, pdfname, orfseqname, protseqname, motif_table){
    
    ### prepare exon information (sequences)
    
    exons_fasta <- read.table(file = exons_fasta, stringsAsFactors = FALSE) # load fasta file
    exons <- exons_fasta$V1[(seq(1, (length(rownames(exons_fasta))-1), 2))] # get first row (exon names, SHOULD BE NUMERIC AND COHERENT WITH THE UNIQ IDs)
    exons_num <- c()
    for (i in exons){ # remove the > & the factor   
        exons_num <- c(exons_num, as.numeric(paste0(strsplit(i, '')[[1]][-1], collapse = '')))
    }
    fastas <- as.character(exons_fasta$V1[(seq(2,length(rownames(exons_fasta)),2))]) # get fastasequences
    
    exons_fasta <- as.data.frame(rbind(exons_num, fastas) ) # merge first and second row together
    colnames(exons_fasta) <- exons_num
    
    ## since I am unable to remove the factors easily, I export and reimport it with stringsAsFactors (hopefully noone else ever reads this)
    
    #write.table(file = "xxx_tmp_pls_delete_xxx", exons_fasta, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    #exons_fasta <- read.table(file = "xxx_tmp_pls_delete_xxx", sep = "\t", stringsAsFactors = FALSE)
    
    # new variants
    
    exons_fasta <- as.data.frame(taRifx::remove.factors(exons_fasta))
    
    
    ### prepare transcript table
    tr.gtf <- read.table(file = conc_tr.gtf, sep = "\t", stringsAsFactors = FALSE) # load concat tr information
    
    # check if it is a reverse gene
    if(unique(tr.gtf$V7) == "-"){ # if strandedness in input is "-"
        reverse_gene <- as.character("yes") # set rev = TRUE
    }else if(unique(tr.gtf$V7)=="+"){
        reverse_gene <- as.character("no")#everything is fine
    }else{print("CAREFUL; DIFFERENT STRANDEDNESS FOR DIFFERENT TRANSCRIPT")}
    
    
    tr_names <- c() 
    for (i in c(1:length(rownames(tr.gtf)))){      # get transcript names
        if (tr.gtf$V3[i]=="transcript"){
            tr_split <- (strsplit(tr.gtf$V9[i], ";"))[[1]][2]   # get tr name
            tr_names <- c(tr_names, gsub(' transcript_id ','',tr_split))
        }
    }
    
    tr <- data.frame(tr_names,c(NA),c(NA),c(NA),c(NA), stringsAsFactors = FALSE) # make new dataframe
    colnames(tr) <- c("V1","V2","V3","V4","V5")
    
    ### merge transcript sequences toegether
    
    for (i in c(1:length(rownames(tr)))){ # for every transcript
        uniqs <- c()
        for (x in (c(1:length(rownames(tr.gtf))))) {# for every possible row (in tr.gtf)
            if(tr.gtf$V3[x]=="exon"){ # for every exon
                tr_split <- gsub(' transcript_id ','',((strsplit(tr.gtf$V9[x],";"))[[1]][2]))  # get belonging transcript    
                if(tr_split == tr$V1[i]){ #if it is the same --> ad uniq identifier
                    uniq_split <- (strsplit(tr.gtf$V9[x], ";"))[[1]][5]
                    #uniqs <- c(uniqs, as.character(gsub(' uniq_id ','', uniq_split)))
                    uniqs <- as.character(paste(uniqs, as.character(gsub(' uniq_id ','', uniq_split)), sep = " "))
                }   
            }                      
        }
        tr$V2[i] <- uniqs       # add exon sequence as input for second coloumn
    }
    
    # add sequence, orf, translation & plots
    orfless <- c()
    for (i in c(1:(length(rownames(tr))))){     # for every transcript
        
        # in the beginning, open pdf
        if(i ==1){
            pdf(file = pdfname)
            par( mfrow = c(2,1))
        }
        
        exons_i <- as.character((tr[i,2]))  # extract the composition of exons
        exons_i <- (strsplit(exons_i, " "))
        exons_i <- as.numeric(exons_i[[1]])
        seq <- c()
        
        for (x in exons_i){                         # for every exon in the transcript
            #print(as.character(exons_fasta[2,x]))        
            seq_x <- as.character(exons_fasta[2,x]) # get sequence
            seq <- c(seq, seq_x)                # add sequence           
        }
        seq <- paste(seq, collapse="")     # collapse
        
        
        if (reverse_gene == "yes"){ # if it is a minus strand
            seq <- DNAString(seq)
            seq <- reverseComplement(seq) # make reverse complement
        }
        
        tr$V3[i] <- tolower(as.character(seq))
        # if no ORF is found
        
        if (is.na(findORFsinSeq(tr$V3[i])[[3]][1])){
            tr$V4[i] <- as.character("NO ORF FOUND")
            tr$V5[i] <- as.character("NO ORF FOUND")
            orfless <- c(orfless,i)
        }else{
            
            orf_dt <- as.data.frame(findORFsinSeq(tr$V3[i]))    # get possible ORFS
            colnames(orf_dt) <- c("start", "stop", "length")    # rename the coloumns
            orf_dt <- orf_dt[order(orf_dt$length, decreasing = TRUE),]  # order by decreasing
            orf_start <- orf_dt$start[1]      # longest is in the first rows
            orf_end <- orf_dt$stop[1]
            tr$V4[i] <- c2s((s2c(tr$V3[i]))[orf_start:orf_end])
            tr$V5[i]<- c2s(seqinr::translate(s2c(tr$V4[i])))
            
        }
        
        if (i == length(rownames(tr))){   # after the last loop
            
            for (n in c(1:(length(rownames(tr))))){   #plot everything
                
                if (n %in% c(orfless)){
                    # don't plot anything
                }else{
                    x_max <- max(nchar(tr$V3))# get longest sequence
                    #plotPotentialStartsAndStopsChangeX(tr$V3[n], xscale = x_max)
                    #title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                    plotORFsinSeqChangeX(tr$V3[n], xscale = x_max)
                    title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                    plotMotifsInSequence(tr$V5[n], motif_table = motif_table)
                    title(main = as.character(tr$V1[n]), font.main = 3, line = 0, outer = FALSE)
                }
            }
            
            
            dev.off()
            output <- c(tr$V4)           # output orfs   
            output <- as.list(output)
            names(output) <- c(as.character(tr$V1))
            write.fasta(sequences = output, names = names(output), file.out = orfseqname)
            
            output <- c(tr$V5)              # output AA sequence
            output <- as.list(output)
            names(output) <- c(as.character(tr$V1))
            
            write.fasta(sequences = output, names = names(output), file.out = protseqname)
        }
        
    }
    
    # add peptide mass
    tr <- cbind(tr,NA)
    colnames(tr)[6]<-"V6"
    
    for(i in c(1:length(rownames(tr)))){
        tr$V6[i]<-CalculatePeptideMass(peptide_sequence = tr$V5[i], unit = "kDa")
    }
    
    # add transcript abundance
    tr <- cbind(tr,NA)
    colnames(tr)[7]<-"V7"
    
    trans_data <- tr.gtf[tr.gtf$V3 == "transcript",]
    
    
    for(i in c(1:length(rownames(tr)))){
        trans_name <- tr$V1[i]
        for(n in allrows(trans_data)){
            tr_tmp <- substr(uncollapse(trans_data$V9[n])[[1]][4], 1, nchar(uncollapse(trans_data$V9[n])[[1]][4])-1)
            if (tr_tmp == trans_name){
                tr$V7[i] <- paste(uncollapse(trans_data$V9[n])[[1]][5:10],collapse=" ")
            }
        }
        
    }
    
    # add in the 8th coloumn: the orf that starts with the signal peptide
    s_pept <- NA
    
    for(i in allrows(motif_table)){
        if(motif_table$V2[i]%in%c("s-pept","S-Pept","signal peptide","spept","Signalpeptid","s pept","s.-pept.","S.-Pept.","S-Peptide","S-Peptid","S.-Peptid","S.-peptide")){
            s_pept<- motif_table$V1[i]
        }
    }
    
    if(is.na(s_pept)){
        print("no signal peptide in motif input")
        # add nothing to 8th coloumkn
        tr <- taRifx::remove.factors(as.data.frame(cbind(tr,"ANALYSIS NOT POSSIBLE")))
        
    }else{ # look @ every ORF --> does it have the S-Pept? --> if yes --> put it into 8th coloumn, if not look for other orfs
        
        col_8 <- c()    
        
        for(i in allrows(tr)){
            print(i)
            #i <- 33
            
            if(grepl(pattern = s_pept, x = tr$V5[i])){ # if s-pept is already in there
                col_8 <- c(col_8, tr$V5[i])
            }else if(tr$V5[i]=="NO ORF FOUND"){
                col_8 <- c(col_8, tr$V5[i])
            }else{
                orf_dt <- as.data.frame(findORFsinSeq(tr$V3[i]))    # get possible ORFS
                colnames(orf_dt) <- c("start", "stop", "length")    # rename the coloumns
                orf_dt <- orf_dt[order(orf_dt$length, decreasing = TRUE),]  # order by decreasing
                orf_dt <-orf_dt[orf_dt$length>=TripleLength(s_pept),] # remove every orf shorter than 3x s-pept
                orf_dt <- orf_dt[-1,] # remove longest (we already used it)
                if(nrow(orf_dt)<1){
                    col_8 <- c(col_8,"-")
                }else{
                    s_pept_found <- FALSE
                    n <- 0
                    while(s_pept_found==FALSE){
                        n <- n + 1
                        orf_start <- orf_dt$start[n]      
                        orf_end <- orf_dt$stop[n]
                        coding_seq <- c2s((s2c(tr$V3[i]))[orf_start:orf_end])
                        pept_seq <- c2s(seqinr::translate(s2c(coding_seq)))
                        if(grepl(pattern = s_pept, x = pept_seq)){
                            col_8 <- c(col_8, pept_seq)
                            tr$V1[i] <- paste(tr$V1[i],"_(mo)",collapse = "", sep  ="")
                            tr$V4[i] <- paste(tr$V4[i],coding_seq, paste=" ", collapse=" ")
                            s_pept_found <- "yes"
                            print("yes")
                        }else if(n == nrow(orf_dt)){ # if the last row is reached and still no orf with the s-pept is found -->
                            col_8 <- c(col_8,"-")
                            s_pept_found <- "no"
                            print("no")
                        }
                    }
                    
                }
            }
            
        }
        
        tr <- taRifx::remove.factors(as.data.frame(cbind(tr,col_8)))
        
        
    }
    
    
    colnames(tr)[8] <- "V8"
    return(tr)
}
# included case of no found ORFs

MergeDuplicateORFs <- function(tr_table){
    
    # tr_table has to be the output of a GetORFAndProteinSeqFrom... - function
    
    tr_orf <- data.frame(tr_table[,c(1,5,7)], stringsAsFactors = FALSE)  # extract only orf and transcript id
    colnames(tr_orf) <- c("tr_names","orf_seqs","quantification")
    # output what happens
    print(paste("Transcripts input: ",length(tr_orf$orf_seqs)," - merging to ",length(unique(tr_orf$orf_seqs))," different transcripts - duplicates: ",(length(tr_orf$orf_seqs)-length(unique(tr_orf$orf_seqs))),".", sep = ""))
    # extract all unique orfs & make new df
    
    orfs_uniq <- unique(tr_orf$orf_seqs)
    
    orf_table <- data.frame(orfs_uniq, c(NA), c(NA),c(NA),c(NA), stringsAsFactors = FALSE)
    colnames(orf_table) <- c("ORF", "Transcript(s)", "frequency", "FPKM_sum", "FPKM_mean")
    
    # find orfs and give tr ids belonging to each orf --> then give number..s
    for (i in c(1:(length(rownames(orf_table))))){ # for every ORF
        transcripts <- c() # clear transcripts string
        orf_i <- orf_table$ORF[i]
        orf_i_df <- tr_orf[tr_orf$orf_seqs==orf_i,] # get dataframe with all trasncript with this ORF
        transcripts <- orf_i_df$tr_names
        orf_table$frequency[i] <- length(transcripts) # add frequency to table
        orf_table$`Transcript(s)`[i]<- paste(transcripts, collapse = " ") # add transcripts to table
        fpkm <- 0
        for(n in transcripts){
            fpkm_tr <- as.numeric(RemoveElements(uncollapse(tr_orf[as.character(tr_orf$tr_names)==n,3])[[1]][4],amount=1))
            fpkm <- fpkm+fpkm_tr
        }
        orf_table$FPKM_sum[i] <- fpkm
        orf_table$FPKM_mean[i] <- orf_table$FPKM_sum[i]/orf_table$frequency[i]
        
    }
    
    # sort orf table by length
    orf_table <- orf_table[order(nchar(orf_table$ORF), decreasing = TRUE),]
    
    return(orf_table)
}
# included FPKM calculation and made input easier

ORFInterpretation <- function(transcripts){
    
    # gene
    gene_name <- transcripts$V1[1]
    while(grepl("_", x = gene_name)){
        gene_name <- RemoveElements(gene_name,amount = 1) 
    }
    
    # orfs
    
    tr_orf <- data.frame(transcripts[,c(1,5,7,8)], stringsAsFactors = FALSE)  # extract only orf and transcript id
    
    tr_orf_new <- c()
    
    for(i in allrows(tr_orf)){
        
        if(grepl(pattern = "(mo)", x = tr_orf$V1[i])){
            line_1 <- tr_orf[i,c(1,2,3)]
            line_1[1] <- RemoveElements(line_1[1],5)
            line_2 <- tr_orf[i, c(1,4,3)]
            line_2[1] <- RemoveElements(line_2[1],5)
            line_2[1] <- paste(line_2[1],"_alt",sep="",paste="")## add name indication!!!!
            names(line_2) <- names(tr_orf_new)
            tr_orf_new <- rbind(tr_orf_new,line_1,line_2)
            
        }else{
            tr_orf_new <- rbind(tr_orf_new,tr_orf[i,c(1,2,3)])
        }
    }
    
    
    
    colnames(tr_orf_new) <- c("tr_names","orf_seqs","quant")
    
    #
    
    # output what happens
    print(paste("Transcripts input: ",length(tr_orf_new$orf_seqs)," - merging to ",length(unique(tr_orf_new$orf_seqs))," different transcripts - duplicates: ",(length(tr_orf_new$orf_seqs)-length(unique(tr_orf_new$orf_seqs))),".", sep = ""))
    # extract all unique orfs & make new df
    
    orfs_uniq <- unique(tr_orf_new$orf_seqs)
    
    orf_table <- data.frame(orfs_uniq, c(NA), c(NA),c(NA),c(NA), stringsAsFactors = FALSE)
    colnames(orf_table) <- c("ORF", "Transcript(s)", "frequency", "FPKM_sum", "FPKM_mean")
    
    # find orfs and give tr ids belonging to each orf --> then give number..s
    for (i in c(1:(length(rownames(orf_table))))){ # for every ORF
        #i <- 1
        orf_transcripts <- c() # clear transcripts string
        orf_i <- orf_table$ORF[i]
        orf_i_df <- tr_orf_new[tr_orf_new$orf_seqs==orf_i,] # get dataframe with all trasncript with this ORF
        orf_transcripts <- orf_i_df$tr_names
        orf_table$frequency[i] <- length(orf_transcripts) # add frequency to table
        orf_table$`Transcript(s)`[i]<- paste(orf_transcripts, collapse = " ") # add transcripts to table
        fpkm <- 0
        for(n in orf_transcripts){
            fpkm_tr <- as.numeric(RemoveElements(uncollapse(tr_orf_new[as.character(tr_orf_new$tr_names)==n,3])[[1]][4],amount=1))
            fpkm <- fpkm+fpkm_tr
        }
        orf_table$FPKM_sum[i] <- fpkm
        orf_table$FPKM_mean[i] <- orf_table$FPKM_sum[i]/orf_table$frequency[i]
        
    }
    
    # sort orf table by length
    orf_table <- orf_table[order(nchar(orf_table$ORF), decreasing = TRUE),]
    
    return(orf_table)
    
    
}
# basically ORF interpretation but now accepts two ORFs for one transcript

MergeDuplicateTranscripts <- function(transcripts, redundant_utr_exons, only_start = TRUE){
    
    # loop over exons and replace redundant start and ends
    
    
    redundant_utr_exons_all <- as.numeric(unlist(uncollapse(redundant_utr_exons)))
    
    for (i in allrows(transcripts)){
        #i <- 2
        comp_i <- uncollapse(transcripts$V2[i])[[1]]
        
        if(as.numeric(comp_i[1]) %in% c(1:100000)){
            # everything fine
        }else{comp_i <- comp_i[-1]} # since sometimes the first is just an empty "" or NA
        
        if(only_start == TRUE){
            n = 1
            ex <- as.numeric(comp_i[n])
            
            if(ex %in% redundant_utr_exons_all){
                for(x in redundant_utr_exons){
                    #x = redundant_utr_exons[5]
                    red_group <- as.numeric(uncollapse(x)[[1]])
                    if(ex %in% red_group){
                        ex <- red_group[1]
                    }
                }
                
                
            }
            
            comp_i[n] <- ex
        }else{
        
        for(n in 1:length(comp_i)){
            #n = 1
            ex <- as.numeric(comp_i[n])
            
            if(ex %in% redundant_utr_exons_all){
                for(x in redundant_utr_exons){
                    #x = redundant_utr_exons[5]
                    red_group <- as.numeric(uncollapse(x)[[1]])
                    if(ex %in% red_group){
                        ex <- red_group[1]
                    }
                }
                
                
            }
            
            comp_i[n] <- ex
        }}    
        
        transcripts$V2[i] <- paste(comp_i,collapse = " ")
        
    }
    
    
    
    #print(paste("Transcripts input: ",length(tr_exons$exon_comp)," - merging to ",length(unique(tr_exons$exon_comp))," different transcripts - duplicates: ",(length(tr_exons$exon_comp)-length(unique(tr_exons$exon_comp))),".", sep = ""))
    trex_uniq <- unique(transcripts$V2)
    
    tr_ex_table <- data.frame(trex_uniq,c(NA),c(NA),c(NA),stringsAsFactors = FALSE)
    colnames(tr_ex_table) <- c("Exon_Composition", "Transcript(s)", "frequency","FPKM_sum")
    # find transcript with same exon composition and give tr_id and number
    for (i in c(1:length(rownames(tr_ex_table)))){ # for every unique Tr_ID
        tr <- c() # clear transcripts string # i = 1
        tr_i <- tr_ex_table$Exon_Composition[i]
        tr_i_df <- transcripts[transcripts$V2==tr_i,] # get dataframe with all trasncript with this composition
        tr <- tr_i_df$V1
        tr_ex_table$frequency[i] <- length(tr) # add frequency to table
        
        fpkm_val <- 0
        for(n in allrows(tr_i_df)){
            # n = 1
            fpkm_val <- fpkm_val + as.numeric(RemoveElements(amount = 1, char = uncollapse(tr_i_df$V7[n])[[1]][4]))
        }
        
        tr_ex_table$FPKM_sum[i] <- fpkm_val
        tr_ex_table$`Transcript(s)`[i]<- paste(tr, collapse = " ") # add transcripts to table
    }
    # sort by exon composition length
    tr_ex_table <- tr_ex_table[order(nchar(tr_ex_table$Exon_Composition), decreasing = TRUE),]
    
    # output
    return(tr_ex_table)
}
# completely new function: give tr_ / transcripts and a string with collapsed redundant utr exons
PlotTranscriptsInExonGraph <- function(exon_data, transcripts, transcript_rows, gene_name, vis.ref.file = FALSE, domains = FALSE, reference = FALSE){
    
    
    ### load data & give new name
    
    exons <- exon_data
    tr <- transcripts
    tr_rows <- transcript_rows
    ref_save <- reference
    
    ref_place <- 0
    suppressWarnings(if(reference == FALSE){
    }else{ # if we have references
        reference <- GetRefTranscripts(ref_table = reference, exons)
        ref_place <- nrow(reference) # make place for them by lowering the number of plotted transcripts
        if(length(tr_rows)>=8){
            tr_rows <- tr_rows[1:(length(tr_rows)-ref_place)]}
    })
    
    
    # wether domains input is given is checked later, but the input should be a dataframe with the first coloumn
    # being the domain name and the second being the exons (in paste(x,y,z,sep=" ") - format)
    
    
    ### parameters of the plot
    
    # check if it is a reverse-stranded gene and multiply the parameters with -1
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else{mult <- as.numeric(1)}
    
    
    ylim_num <- 1300 # highest point of plot
    
    ref <- c(50+900,100+900) # first line
    addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
    max_left <- 100 * mult # leftmost coordinate
    max_right <- 10100 * mult# rightmost coordinate
    
    cols_coverage <- colorRampPalette(colors = c("black","#E65C5C", "darkgreen")) # colors - first is for lowly abundant exons, last is for highly abundant exons
    # cols_coverage(6)[1:6] 
    col_ref <- "darkblue"
    ##E65C5C = red
    ### exon preprocessing
    
    # since changing the exon coordinates each time this function is called for the same dataset takes long time (2-3s), I implemented
    # the possibility of saving the coordinate-dataset into .env
    # first, the function checks if there is a dataset already existing for this exon set (based on name)
    # then, it checks wether the length of the rows is the same --> if true --> use this dataset and skip directly to the plotting
    # if vis.ref == NA (default) or there is no, or if the rows don't match --> do the changing and save a new vis.ref
    
    
    
    if (vis.ref.file == FALSE){ # if no vis.ref is given --> set vis.ref.yes to no
        vis.ref.yes <- "no"
    }else{vis.ref.yes <- "yes" # else set to yes, but
    
    # the next expression is complicated, so ill explain it here in-depth:
    # I want to read.table vis.ref.file and make the table vis.ref... but if vis.ref.file does not exist yet, (but shall exist after this function worked once), the function
    # gives me an error. Thats why I made that statement into the tryCatch function, which gives makes the following: instead of an error message, it checks wether there is an "old" vis.ref in the environment
    # and if yes, it deletes that, so that no vis.ref is anywhere
    tryCatch(expr = vis.ref <-read.table(file = vis.ref.file, sep = "\t", stringsAsFactors = FALSE),
             warning=function(error_message){if(exists("vis.ref")){rm(vis.ref, envir = .GlobalEnv)}} )
    if(exists("vis.ref")){# check if it exists
        if(nrow(exons) == nrow(vis.ref)){# check wether it has the same number of rows
            exons <- as.data.frame(vis.ref) # overwrite exons with the vis.ref
        }else{vis.ref.yes <- "no"} # if not --> vis.ref.yes <- no                     
    }else{vis.ref.yes <- "no"} 
    
    }
    
    if(vis.ref.yes == "no"){ # if no vis.ref is given do the analysis
        
        # find first exon beginning and substract this value
        
        min_start <- min(exons$V4)
        exons$V4 <- exons$V4 - min_start
        exons$V5 <- exons$V5 - min_start
        
        # find overlap groups
        
        FindFirstOverlapGroup <- function(exon_df, rows){
            df_ex <- as.data.frame(exon_df)
            for (i in c(rows)){   # for every exon
                if (i == min(rows)){    # if it is the first exon of the input
                    group <- c(i)
                    range_a <- c(df_ex$V4[i]:df_ex$V5[i])   # define first range
                }else{      # for every other exon
                    range_tmp <- c(df_ex$V4[i]:df_ex$V5[i]) # define its own range (_tmp)    
                    is_part_of_a <- 0         # set this parameter to 0
                    for (n in range_tmp){   # for ever value in range_tmp
                        if(n%in%range_a){   # check if it is in range_a
                            is_part_of_a <- is_part_of_a+1  # and if this is the case, increment is_part_of by one 
                        }
                    }
                    if (is_part_of_a >= 1){   # if is_part_of >= 1, this means, that range_temp is part of range_a
                        range_a <- unique(c(range_a, range_tmp)) # merge them 
                        group <- c(group,i)   # and add to group
                    } else{
                        # not_in_group <- c(not_in_group,i)
                    }
                }
            }
            #print(group)
            return(group)
        }
        
        # group <- 18:21
        firsts <- c(0) # make sure to delete the zero later --> zero gets removed somehow automatically, yeah!
        rowsmax <- length(rownames(exons))
        longest <- c(0)
        #i <- 2
        for (i in c(1:rowsmax)){ # for every row, (theoretically, it ends before, but it COULD be the case, that there are rowmax different exons with no overlap overall)
            if (max(firsts)<length(rownames(exons))){#if the highest/last exon added to first is still lower than the last exon --> do the analysis
                if (i ==1){
                    group <- FindFirstOverlapGroup(exons,rows = (c(1:rowsmax)))
                    firsts <- group[length(FindFirstOverlapGroup(exons,rows = (c(1:rowsmax))))] # get the last of the first overlap group
                    group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                    longest <- rownames(group_tmp[which.max(group_tmp$V5),])[1] # if there are multiple --> take the first, it doesnt matter                
                    
                }else{
                    group <- FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax)))
                    firsts <- c(firsts,group[length(FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax))))]) # get the last exon of the next overlap group
                    group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                    longest <- as.numeric(c(longest,rownames(group_tmp[which.max(group_tmp$V5),])[1]))                
                    
                }
            }}
        # its okay if there is an error message
        
        #firsts # these are all the last exons of the overlap group
        firsts <- firsts[-length(firsts)] # remove the last one
        seconds <- firsts+1 # these are all exons whose are the first exons of the following overlap group
        longest <- as.numeric(longest[-length(longest)])
        
        exon_dist_rows <- as.data.frame(cbind(firsts,seconds,longest)) # merge to dataframe
        
        subtrExonDist <- function(exon_table,before,after,longest){ # make function
            subtr <- exon_table$V4[after]-exon_table$V5[longest] - 100   # get intron length - 100
            exon_table$V4[c(after:length(rownames(exon_table)))] <- exon_table$V4[c(after:length(rownames(exon_table)))]-subtr # subtract this value
            exon_table$V5[c(after:length(rownames(exon_table)))] <- exon_table$V5[c(after:length(rownames(exon_table)))]-subtr # subtract this value 
            return(exon_table) # return the new exon table
        }
        
        #i <- 1
        for (i in c(1:length(rownames(exon_dist_rows)))){ # for every intron beginning
            bf <- exon_dist_rows$firsts[i]  # get start number
            af <- exon_dist_rows$seconds[i] # get end number
            lon <- exon_dist_rows$longest[i]
            af_minus_long <- as.numeric(exons$V4[af])-as.numeric(exons$V5[lon]) # check distance/intron length
            if (as.numeric(af_minus_long) >= 100){ # if bigger or equal 100
                exons <- subtrExonDist(exon_table = exons, before = bf, after = af, longest = lon) # transform distance to 100
            }
            
        }
        
        
        # normalize for possible place: 10000 Units
        
        max_end <- max(exons$V5) # get highest
        
        fact <- 10000/max_end # get multiplicator
        exons$V4 <- fact*exons$V4 # multiplicate
        exons$V5 <- fact*exons$V5
        exons$V4 <- round(digits = 1,x = exons$V4) # round the results
        exons$V5 <- round(digits = 1, x = exons$V5)
        
        if(vis.ref.file == FALSE){ # do nothing
        }else{ # make new vis.ref.file and save it 
            write.table(x = exons, file = vis.ref.file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
        }
        
    } # the closing from if vis.rev.yes == no
    
    exons$V4 <- exons$V4*mult
    exons$V5 <- exons$V5*mult
    
    ### plotting
    #alex11()
    
    if(mult==as.numeric(-1)){
        additional_place <- 200
    }else{additional_place <- 0}
    
    plot_x1 <- (additional_place+10200) * mult
    plot_x2 <- -300 * mult
    plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
    #rect(xleft = 1, xright = 2, ybottom = 0, ytop = 1300, border = "NA", col = "gray")  # plot "y-axis"
    
    # introns
    rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
    
    # plot every possible exon in the highest line
    
    for (i in c(1:(length(rownames(exons))))){
        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
    }
    
    # check if domains are specified and if yes, add them below the plot
    
    suppressWarnings(if(domains == FALSE){
        text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
    }else{
        
        
        dom <- domains
        
        for(i in allrows(dom)){
            #i <- 5
            ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
            
            
            firstexon <- floor(ex_seq[1]) # get integer exon
            perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
            left_coord <- exons$V4[firstexon]+((exons$V5[firstexon]-exons$V4[firstexon])*perc_f_exon)
            
            
            lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
            perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
            if(ex_seq[length(ex_seq)]%%1 > 0){
                right_coord <- exons$V4[lastexon]+((exons$V5[lastexon]-exons$V4[lastexon])*perc_l_exon)
            }else{right_coord <- exons$V5[lastexon]}
            
            ex_seq[1] <- firstexon
            ex_seq[length(ex_seq)] <- lastexon
            if(length(unique(ex_seq))==1){
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }else{for(n in ex_seq){
                if(n==firstexon){
                    rect(xleft = max_left+left_coord, xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else if(n==lastexon){
                    rect(xleft = max_left+exons$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else{rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }
                
                
                #rect(xleft = max_left+left_coord, xright =max_left+left_coord+10, ybot = ref[1]-1000, ytop = ref[2]-1000, col = "grey", border = TRUE)# remove the black line after
                #rect(xleft = max_left+right_coord-10, xright =max_left+right_coord, ybot = ref[1]-1000, ytop = ref[2]-1000, col = "grey", border = TRUE)# remove the black line after
                
                mean_dom <- max_left + (left_coord+right_coord)/2
                #text(labels = (paste(dom[i,1],sep="")), x = mean_dom+(1/2*(75*nchar(c2s(dom$domains[1])))), y =ref[2]-1075, cex = 0.8,  srt=45)
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                #rect(xleft = mean_dom, xright = mean_dom, ybot = ref[1]-1000, ytop = ref[1], col = "blue", border = TRUE)# remove the black line after
                
            }}
            
            #if(length(ex_seq)>=2){   
            rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
            #}
            
            
        }
    })
    
    # plot transcripts
    
    # if the strand is "-" --> turn around / reverse sequence and orf (for orf detection et cetera)
    #strand <- exons$V7[1]
    #if(strand == "-"){
    
    #    for (i in allrows(tr)){
    #i <- 1
    #        seq <- tr$V3[i]
    #        orf <- tr$V4[i]
    #        tr$V3[i] <- paste(rev(substring(seq,1:nchar(seq),1:nchar(seq))),collapse="")
    #        tr$V4[i] <- paste(rev(substring(orf,1:nchar(orf),1:nchar(orf))),collapse="") 
    
    #    }
    #}
    
    
    # n <- 3
    for (n in tr_rows){  # for chosen transcript rows
        ex_seq <- uncollapse(tr$V2[n])[[1]][-1]  #get exon seuqnece
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        cov_i <- uncollapse(tr$V7[n])[[1]][4]    # uncollapse and extract FPKM for this transcript
        cov_i <- as.numeric(substr(cov_i, 1, nchar(cov_i)-1))
        
        if (cov_i < 0.3){ # get colour depending on coverage
            col_cov <- cols_coverage(6)[1]
        } else if (cov_i < 1){
            col_cov <- cols_coverage(6)[2]
        } else if (cov_i < 4){
            col_cov <- cols_coverage(6)[3]
        } else if (cov_i < 6){
            col_cov <- cols_coverage(6)[4]
        } else if (cov_i < 10){
            col_cov <- cols_coverage(6)[5]
        } else { col_cov <- cols_coverage(6)[6]}
        
        
        if (n == tr_rows[1]){   # get x as the linenumber of this transcript
            x <- 1 + ref_place
        }else if( n == tr_rows[2]){
            x <- 2 + ref_place
        }else if( n == tr_rows[3]){
            x <- 3 + ref_place
        }else if (n == tr_rows[4]){
            x <- 4 + ref_place
        }else if(n == tr_rows[5]){
            x <- 5 + ref_place
        }else if(n == tr_rows[6]){
            x <- 6 + ref_place
        }else if(n == tr_rows[7]){
            x <- 7 + ref_place
        }else if(n == tr_rows[8]){
            x <- 8 + ref_place
        }else {x <- 9 + ref_place}
        
        # get exons which are in the ORF and in the UTR
        
        if(grepl(x = tr$V1[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
            orf <- as.character(tr$V4[n]) # proceed as usual
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            if(mult==as.numeric(-1)){
                text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            }else{
                text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            }
            
            
        }else{ # if there is one
            orf <- as.character(uncollapse(tr$V4[n])[[1]][1]) # give him the first (longest)
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            if(mult==as.numeric(-1)){
                text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            }else{
                text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            }
            
            ################# # give the second orf!!!
            
            
            
            
            orf <- as.character(uncollapse(tr$V4[n])[[1]][2]) # give him the first (longest)
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            
            
            
            
            
            
        }
        
        
    }
    
    
    # add reference
    x <- 0
    suppressWarnings(if(reference == FALSE){
    }else{
        col_cov <- "deepskyblue"
        for(i in allrows(reference)){ # i <- 1
            ref_tr_name <- reference$V1[i]
            x <- x+1 # next line for plotting the reference
            ex_seq <- uncollapse(reference$V2[i])[[1]] #get exon seuqnece
            if("?" %in% ex_seq){
                # print("reference with unknown exon skipped")
                # else (it is ? ) --> look at exon_table (not transformed)
                # check wether start is found anywhere > save value from transformed exons
                # if not: loop through every exon and check wether it lies in an exon
                # if yes: get proportion and multiply with value from transformed exon > save
                # if not: check wether it lies between introns (careful with first and last line)
                # if yes: get proportion and multiply with value from transformed exon > save
                
                # same for end
                # plot as normal
                
                # get all exons of this transcript
                ref_tr <- c()
                for(y in allrows(ref_save)){ # look @ every row y<-3
                    if(ref_save$V3[y]=="exon"){ # only look @ exons
                        if(ref_tr_name == (RemoveElements(char = (uncollapse(ref_save$V9[y])[[1]][4]), amount = 1, start = "last"))){
                            ref_tr <- as.data.frame(rbind(ref_tr,ref_save[y,]))
                        }
                    }
                }
                
                
                #x <- 2
                position <- 0    
                for(n in ex_seq){ # loop over ex_seq             # n <- 39
                    position <- position + 1
                    suppressWarnings(n <- as.numeric(n))
                    if(is.na(n)){ # if n was not a numeric value (? becomes NA because as.numeric)
                        start_n <- ref_tr$V4[position] # get start and stop value for this exon 
                        stop_n <- ref_tr$V5[position]
                        
                        ### look at start value
                        exonic_start <- "nay"
                        inbetween_exon <- "nay"
                        for (exons_orig in allrows(exon_data)){ # loop over every exon
                            if(exon_data$V4[exons_orig]==start_n){ # if one exons has the same start
                                start_ex <- exons_orig # define this as the start exon
                                exonic_start <- "yay"
                                
                            }
                        } 
                        if (exonic_start == "yay"){
                            start_n <- exons$V4[start_ex] # assign the start point from the table of the transformed exons
                        }else{ # the start location is not known :( --> need to find the region to approximate
                            for (exons_orig in allrows(exon_data)){ # check first if it is between 
                                if(exon_data$V4[exons_orig]<start_n & exon_data$V5[exons_orig]>start_n){ # if it lies inbetween an exon 
                                    perc <- (start_n-exon_data$V4[exons_orig])/(exon_data$V5[exons_orig]-exon_data$V4[exons_orig])
                                    start_n <- exons$V4[exons_orig]+((exons$V5[exons_orig]-exons$V4[exons_orig])*perc)
                                    inbetween_exon <- "yay"
                                }
                                
                            }
                            
                            if (inbetween_exon == "nay"){
                                for (exons_orig in allrows(exon_data)){ # check first if it is between
                                    if(exons_orig == max(allrows(exon_data))){
                                        
                                    }else{
                                        if(exon_data$V5[exons_orig]<start_n & exon_data$V4[exons_orig+1]>start_n){ # if it lies inbetween an intron 
                                            perc <- (start_n-exon_data$V5[exons_orig])/(exon_data$V4[exons_orig+1]-exon_data$V5[exons_orig])
                                            start_n <- exons$V5[exons_orig]+((exons$V4[exons_orig+1]-exons$V5[exons_orig])*perc)
                                        }
                                    }  
                                } 
                            } 
                            
                            
                        }
                        
                        ###
                        ### look at end value
                        exonic_end <- "nay"
                        inbetween_exon_end <- "nay"
                        for (exons_orig in allrows(exon_data)){ # loop over every exon
                            if(exon_data$V5[exons_orig]==stop_n){ # if one exons has the same start
                                stop_ex <- exons_orig # define this as the start exon
                                exonic_end <- "yay"
                                
                            }
                        } 
                        if (exonic_end == "yay"){
                            stop_n <- exons$V5[stop_ex] # assign the start point from the table of the transformed exons
                        }else{ # the start location is not known :( --> need to find the region to approximate
                            for (exons_orig in allrows(exon_data)){ # check first if it is between 
                                if(exon_data$V4[exons_orig]<stop_n & exon_data$V5[exons_orig]>stop_n){ # if it lies inbetween an exon 
                                    perc <- (stop_n-exon_data$V4[exons_orig])/(exon_data$V5[exons_orig]-exon_data$V4[exons_orig])
                                    stop_n <- exons$V4[exons_orig]+((exons$V5[exons_orig]-exons$V4[exons_orig])*perc)
                                    inbetween_exon_end <- "yay"
                                }
                                
                            }
                            
                            if (inbetween_exon_end == "nay"){
                                for (exons_orig in allrows(exon_data)){ # check first if it is between
                                    if(exons_orig == max(allrows(exon_data))){
                                        
                                    }else{
                                        if(exon_data$V5[exons_orig]<stop_n & exon_data$V4[exons_orig+1]>stop_n){ # if it lies inbetween an intron 
                                            perc <- (stop_n-exon_data$V5[exons_orig])/(exon_data$V4[exons_orig+1]-exon_data$V5[exons_orig])
                                            stop_n <- exons$V5[exons_orig]+((exons$V4[exons_orig+1]-exons$V5[exons_orig])*perc)
                                        }
                                    }
                                } 
                            }
                            
                            
                        }
                        
                        
                        if (position == 1){ # same for first one
                            #max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                            #max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        rect(xleft = max_left+start_n, xright = max_left+stop_n,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                        
                    }else{# if n was a numeric value
                        
                        if (position == 1){ # same for first one
                            max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                        #rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "black",density=15, angle=90, border = FALSE)
                        
                        text(x = max_left+exons$V4[n]+20, y = ref[1]+addtoref[x]-20, labels = n, cex = 0.65)
                        
                        
                        
                    }
                    
                }
            }else{ # if there is no ? in the exon sequence
                for (n in ex_seq){  # loop over exon sequence
                    n <- as.numeric(n)
                    
                    if (n == as.numeric(ex_seq[1])){ # same for first one
                        max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                    #rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "black",density=15, angle=90, border = FALSE)
                    
                    text(x = max_left+exons$V4[n]+20, y = ref[1]+addtoref[x]-20, labels = n, cex = 0.65)
                    
                    
                }
                text(x = max_left - 450*mult, y =  ref[1]+addtoref[x]+25, labels = reference$V1[i]) 
            }
        }
    }
    
    
    )
    
    # add exonic composition of the genetic region 
    
    ex_sorted <- as.data.frame(exons)
    
    addtorefexon <-c(10,27,44,61,78,95,112,129)
    #x <- 7
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                end <- abs(ex_sorted$V5[i])
            }else if((abs(ex_sorted$V4[i])-end)>=15){
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- abs(ex_sorted$V5[i]) 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    
    
    # add legend
    if(mult==as.numeric(-1)){
        xlegend <- 550*mult}else{xlegend <- 7550}
    
    legend(legend = c("< 0.3 FPKM","< 1 FPKM","< 4 FPKM","< 6 FPKM","< 10 FPKM", "> 10 FPKM"), x = xlegend, y = ylim_num+120, fill=c(cols_coverage(6)[1:6]))#, title = "colour legend") 
    
    # add descriptions
    
    if(mult==as.numeric(-1)){
        text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
        
    }else{
        text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        
    }
    
    # add absolute location
    
    chr <- exon_data$V1[1]
    if(mult==as.numeric(-1)){
        gene_start <- min(exon_data$V4)
        gene_end <- max(exon_data$V5)
        gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
        text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
        
        
        
    }else{
        gene_start <- min(exon_data$V4)
        gene_end <- max(exon_data$V5)
        gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
        text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
        
    }
    
    
    
    # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
    
    text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
    
    
    
}
# if two orfs exist, plots two orfs
PlotTranscriptsInExonGraph <- function(exon_data, transcripts, transcript_rows, gene_name, vis.ref.file = FALSE, domains = FALSE, reference = FALSE){
    
    
    ### load data & give new name
    
    exons <- exon_data
    tr <- transcripts
    tr_rows <- transcript_rows
    ref_save <- reference
    
    ref_place <- 0
    suppressWarnings(if(reference == FALSE){
    }else{ # if we have references
        reference <- GetRefTranscripts(ref_table = reference, exons)
        ref_place <- nrow(reference) # make place for them by lowering the number of plotted transcripts
        if(length(tr_rows)>=8){
            tr_rows <- tr_rows[1:(length(tr_rows)-ref_place)]}
    })
    
    
    # wether domains input is given is checked later, but the input should be a dataframe with the first coloumn
    # being the domain name and the second being the exons (in paste(x,y,z,sep=" ") - format)
    
    
    ### parameters of the plot
    
    # check if it is a reverse-stranded gene and multiply the parameters with -1
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else{mult <- as.numeric(1)}
    
    
    ylim_num <- 1300 # highest point of plot
    
    ref <- c(50+900,100+900) # first line
    addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
    max_left <- 100 * mult # leftmost coordinate
    max_right <- 10100 * mult# rightmost coordinate
    
    cols_coverage <- colorRampPalette(colors = c("black","#E65C5C", "darkgreen")) # colors - first is for lowly abundant exons, last is for highly abundant exons
    # cols_coverage(6)[1:6] 
    col_ref <- "darkblue"
    ##E65C5C = red
    ### exon preprocessing
    
    # since changing the exon coordinates each time this function is called for the same dataset takes long time (2-3s), I implemented
    # the possibility of saving the coordinate-dataset into .env
    # first, the function checks if there is a dataset already existing for this exon set (based on name)
    # then, it checks wether the length of the rows is the same --> if true --> use this dataset and skip directly to the plotting
    # if vis.ref == NA (default) or there is no, or if the rows don't match --> do the changing and save a new vis.ref
    
    
    
    if (vis.ref.file == FALSE){ # if no vis.ref is given --> set vis.ref.yes to no
        vis.ref.yes <- "no"
    }else{vis.ref.yes <- "yes" # else set to yes, but
    
    # the next expression is complicated, so ill explain it here in-depth:
    # I want to read.table vis.ref.file and make the table vis.ref... but if vis.ref.file does not exist yet, (but shall exist after this function worked once), the function
    # gives me an error. Thats why I made that statement into the tryCatch function, which gives makes the following: instead of an error message, it checks wether there is an "old" vis.ref in the environment
    # and if yes, it deletes that, so that no vis.ref is anywhere
    tryCatch(expr = vis.ref <-read.table(file = vis.ref.file, sep = "\t", stringsAsFactors = FALSE),
             warning=function(error_message){if(exists("vis.ref")){rm(vis.ref, envir = .GlobalEnv)}} )
    if(exists("vis.ref")){# check if it exists
        if(nrow(exons) == nrow(vis.ref)){# check wether it has the same number of rows
            exons <- as.data.frame(vis.ref) # overwrite exons with the vis.ref
        }else{vis.ref.yes <- "no"} # if not --> vis.ref.yes <- no                     
    }else{vis.ref.yes <- "no"} 
    
    }
    
    if(vis.ref.yes == "no"){ # if no vis.ref is given do the analysis
        
        # find first exon beginning and substract this value
        
        min_start <- min(exons$V4)
        exons$V4 <- exons$V4 - min_start
        exons$V5 <- exons$V5 - min_start
        
        # find overlap groups
        
        FindFirstOverlapGroup <- function(exon_df, rows){
            df_ex <- as.data.frame(exon_df)
            for (i in c(rows)){   # for every exon
                if (i == min(rows)){    # if it is the first exon of the input
                    group <- c(i)
                    range_a <- c(df_ex$V4[i]:df_ex$V5[i])   # define first range
                }else{      # for every other exon
                    range_tmp <- c(df_ex$V4[i]:df_ex$V5[i]) # define its own range (_tmp)    
                    is_part_of_a <- 0         # set this parameter to 0
                    for (n in range_tmp){   # for ever value in range_tmp
                        if(n%in%range_a){   # check if it is in range_a
                            is_part_of_a <- is_part_of_a+1  # and if this is the case, increment is_part_of by one 
                        }
                    }
                    if (is_part_of_a >= 1){   # if is_part_of >= 1, this means, that range_temp is part of range_a
                        range_a <- unique(c(range_a, range_tmp)) # merge them 
                        group <- c(group,i)   # and add to group
                    } else{
                        # not_in_group <- c(not_in_group,i)
                    }
                }
            }
            #print(group)
            return(group)
        }
        
        # group <- 18:21
        firsts <- c(0) # make sure to delete the zero later --> zero gets removed somehow automatically, yeah!
        rowsmax <- length(rownames(exons))
        longest <- c(0)
        #i <- 2
        for (i in c(1:rowsmax)){ # for every row, (theoretically, it ends before, but it COULD be the case, that there are rowmax different exons with no overlap overall)
            if (max(firsts)<length(rownames(exons))){#if the highest/last exon added to first is still lower than the last exon --> do the analysis
                if (i ==1){
                    group <- FindFirstOverlapGroup(exons,rows = (c(1:rowsmax)))
                    firsts <- group[length(FindFirstOverlapGroup(exons,rows = (c(1:rowsmax))))] # get the last of the first overlap group
                    group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                    longest <- rownames(group_tmp[which.max(group_tmp$V5),])[1] # if there are multiple --> take the first, it doesnt matter                
                    
                }else{
                    group <- FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax)))
                    firsts <- c(firsts,group[length(FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax))))]) # get the last exon of the next overlap group
                    group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                    longest <- as.numeric(c(longest,rownames(group_tmp[which.max(group_tmp$V5),])[1]))                
                    
                }
            }}
        # its okay if there is an error message
        
        #firsts # these are all the last exons of the overlap group
        firsts <- firsts[-length(firsts)] # remove the last one
        seconds <- firsts+1 # these are all exons whose are the first exons of the following overlap group
        longest <- as.numeric(longest[-length(longest)])
        
        exon_dist_rows <- as.data.frame(cbind(firsts,seconds,longest)) # merge to dataframe
        
        subtrExonDist <- function(exon_table,before,after,longest){ # make function
            subtr <- exon_table$V4[after]-exon_table$V5[longest] - 100   # get intron length - 100
            exon_table$V4[c(after:length(rownames(exon_table)))] <- exon_table$V4[c(after:length(rownames(exon_table)))]-subtr # subtract this value
            exon_table$V5[c(after:length(rownames(exon_table)))] <- exon_table$V5[c(after:length(rownames(exon_table)))]-subtr # subtract this value 
            return(exon_table) # return the new exon table
        }
        
        #i <- 1
        for (i in c(1:length(rownames(exon_dist_rows)))){ # for every intron beginning
            bf <- exon_dist_rows$firsts[i]  # get start number
            af <- exon_dist_rows$seconds[i] # get end number
            lon <- exon_dist_rows$longest[i]
            af_minus_long <- as.numeric(exons$V4[af])-as.numeric(exons$V5[lon]) # check distance/intron length
            if (as.numeric(af_minus_long) >= 100){ # if bigger or equal 100
                exons <- subtrExonDist(exon_table = exons, before = bf, after = af, longest = lon) # transform distance to 100
            }
            
        }
        
        
        # normalize for possible place: 10000 Units
        
        max_end <- max(exons$V5) # get highest
        
        fact <- 10000/max_end # get multiplicator
        exons$V4 <- fact*exons$V4 # multiplicate
        exons$V5 <- fact*exons$V5
        exons$V4 <- round(digits = 1,x = exons$V4) # round the results
        exons$V5 <- round(digits = 1, x = exons$V5)
        
        if(vis.ref.file == FALSE){ # do nothing
        }else{ # make new vis.ref.file and save it 
            write.table(x = exons, file = vis.ref.file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
        }
        
    } # the closing from if vis.rev.yes == no
    
    exons$V4 <- exons$V4*mult
    exons$V5 <- exons$V5*mult
    
    ### plotting
    #alex11()
    
    if(mult==as.numeric(-1)){
        additional_place <- 200
    }else{additional_place <- 0}
    
    plot_x1 <- (additional_place+10200) * mult
    plot_x2 <- -300 * mult
    plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
    #rect(xleft = 1, xright = 2, ybottom = 0, ytop = 1300, border = "NA", col = "gray")  # plot "y-axis"
    
    # introns
    rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
    
    # plot every possible exon in the highest line
    
    for (i in c(1:(length(rownames(exons))))){
        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
    }
    
    # check if domains are specified and if yes, add them below the plot
    
    suppressWarnings(if(domains == FALSE){
        text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
    }else{
        
        
        dom <- domains
        
        for(i in allrows(dom)){
            #i <- 5
            ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
            
            
            firstexon <- floor(ex_seq[1]) # get integer exon
            perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
            left_coord <- exons$V4[firstexon]+((exons$V5[firstexon]-exons$V4[firstexon])*perc_f_exon)
            
            
            lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
            perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
            if(ex_seq[length(ex_seq)]%%1 > 0){
                right_coord <- exons$V4[lastexon]+((exons$V5[lastexon]-exons$V4[lastexon])*perc_l_exon)
            }else{right_coord <- exons$V5[lastexon]}
            
            ex_seq[1] <- firstexon
            ex_seq[length(ex_seq)] <- lastexon
            if(length(unique(ex_seq))==1){
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }else{for(n in ex_seq){
                if(n==firstexon){
                    rect(xleft = max_left+left_coord, xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else if(n==lastexon){
                    rect(xleft = max_left+exons$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else{rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }
                
                
                #rect(xleft = max_left+left_coord, xright =max_left+left_coord+10, ybot = ref[1]-1000, ytop = ref[2]-1000, col = "grey", border = TRUE)# remove the black line after
                #rect(xleft = max_left+right_coord-10, xright =max_left+right_coord, ybot = ref[1]-1000, ytop = ref[2]-1000, col = "grey", border = TRUE)# remove the black line after
                
                mean_dom <- max_left + (left_coord+right_coord)/2
                #text(labels = (paste(dom[i,1],sep="")), x = mean_dom+(1/2*(75*nchar(c2s(dom$domains[1])))), y =ref[2]-1075, cex = 0.8,  srt=45)
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                #rect(xleft = mean_dom, xright = mean_dom, ybot = ref[1]-1000, ytop = ref[1], col = "blue", border = TRUE)# remove the black line after
                
            }}
            
            #if(length(ex_seq)>=2){   
            rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
            #}
            
            
        }
    })
    
    # plot transcripts
    
    # if the strand is "-" --> turn around / reverse sequence and orf (for orf detection et cetera)
    #strand <- exons$V7[1]
    #if(strand == "-"){
    
    #    for (i in allrows(tr)){
    #i <- 1
    #        seq <- tr$V3[i]
    #        orf <- tr$V4[i]
    #        tr$V3[i] <- paste(rev(substring(seq,1:nchar(seq),1:nchar(seq))),collapse="")
    #        tr$V4[i] <- paste(rev(substring(orf,1:nchar(orf),1:nchar(orf))),collapse="") 
    
    #    }
    #}
    
    
    # n <- tr_rows[1]
    for (n in tr_rows){  # for chosen transcript rows
        ex_seq <- uncollapse(tr$V2[n])[[1]][-1]  #get exon seuqnece
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        cov_i <- uncollapse(tr$V7[n])[[1]][4]    # uncollapse and extract FPKM for this transcript
        cov_i <- as.numeric(substr(cov_i, 1, nchar(cov_i)-1))
        
        if (cov_i < 0.3){ # get colour depending on coverage
            col_cov <- cols_coverage(6)[1]
        } else if (cov_i < 1){
            col_cov <- cols_coverage(6)[2]
        } else if (cov_i < 4){
            col_cov <- cols_coverage(6)[3]
        } else if (cov_i < 6){
            col_cov <- cols_coverage(6)[4]
        } else if (cov_i < 10){
            col_cov <- cols_coverage(6)[5]
        } else { col_cov <- cols_coverage(6)[6]}
        
        
        if (n == tr_rows[1]){   # get x as the linenumber of this transcript
            x <- 1 + ref_place
        }else if( n == tr_rows[2]){
            x <- 2 + ref_place
        }else if( n == tr_rows[3]){
            x <- 3 + ref_place
        }else if (n == tr_rows[4]){
            x <- 4 + ref_place
        }else if(n == tr_rows[5]){
            x <- 5 + ref_place
        }else if(n == tr_rows[6]){
            x <- 6 + ref_place
        }else if(n == tr_rows[7]){
            x <- 7 + ref_place
        }else if(n == tr_rows[8]){
            x <- 8 + ref_place
        }else {x <- 9 + ref_place}
        
        # get exons which are in the ORF and in the UTR
        
        if(grepl(x = tr$V1[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
            orf <- as.character(tr$V4[n]) # proceed as usual
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            if(mult==as.numeric(-1)){
                text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            }else{
                text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            }
            
            
        }else{ # if there is one
            orf <- as.character(uncollapse(tr$V4[n])[[1]][1]) # give him the first (longest)
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            if(mult==as.numeric(-1)){
                text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            }else{
                text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            }
            
            ################# # give the second orf!!!
            
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            
            orf <- as.character(uncollapse(tr$V4[n])[[1]][2]) # give him the first (longest)
            sequence <- as.character(tr$V3[n])    
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            
            #ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            # if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
            #    ex_seq <- ex_seq[length(ex_seq):1]
            #}
            
            # i <- 73
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){ # if the orf lies in only one exon
                    rect(xleft = max_left+orf_start_exon, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)            
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov, border = FALSE)
                    }
                }
                
            }
            
            if(mult==as.numeric(-1)){
                text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            }else{
                text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            }
            
            
            
            
            ################ end of insertion
            
            
            
            
        }
        
        
    }
    
    
    # add reference
    x <- 0
    suppressWarnings(if(reference == FALSE){
    }else{
        col_cov <- "deepskyblue"
        for(i in allrows(reference)){ # i <- 1
            ref_tr_name <- reference$V1[i]
            x <- x+1 # next line for plotting the reference
            ex_seq <- uncollapse(reference$V2[i])[[1]] #get exon seuqnece
            if("?" %in% ex_seq){
                # print("reference with unknown exon skipped")
                # else (it is ? ) --> look at exon_table (not transformed)
                # check wether start is found anywhere > save value from transformed exons
                # if not: loop through every exon and check wether it lies in an exon
                # if yes: get proportion and multiply with value from transformed exon > save
                # if not: check wether it lies between introns (careful with first and last line)
                # if yes: get proportion and multiply with value from transformed exon > save
                
                # same for end
                # plot as normal
                
                # get all exons of this transcript
                ref_tr <- c()
                for(y in allrows(ref_save)){ # look @ every row y<-3
                    if(ref_save$V3[y]=="exon"){ # only look @ exons
                        if(ref_tr_name == (RemoveElements(char = (uncollapse(ref_save$V9[y])[[1]][4]), amount = 1, start = "last"))){
                            ref_tr <- as.data.frame(rbind(ref_tr,ref_save[y,]))
                        }
                    }
                }
                
                
                #x <- 2
                position <- 0    
                for(n in ex_seq){ # loop over ex_seq             # n <- 39
                    position <- position + 1
                    suppressWarnings(n <- as.numeric(n))
                    if(is.na(n)){ # if n was not a numeric value (? becomes NA because as.numeric)
                        start_n <- ref_tr$V4[position] # get start and stop value for this exon 
                        stop_n <- ref_tr$V5[position]
                        
                        ### look at start value
                        exonic_start <- "nay"
                        inbetween_exon <- "nay"
                        for (exons_orig in allrows(exon_data)){ # loop over every exon
                            if(exon_data$V4[exons_orig]==start_n){ # if one exons has the same start
                                start_ex <- exons_orig # define this as the start exon
                                exonic_start <- "yay"
                                
                            }
                        } 
                        if (exonic_start == "yay"){
                            start_n <- exons$V4[start_ex] # assign the start point from the table of the transformed exons
                        }else{ # the start location is not known :( --> need to find the region to approximate
                            for (exons_orig in allrows(exon_data)){ # check first if it is between 
                                if(exon_data$V4[exons_orig]<start_n & exon_data$V5[exons_orig]>start_n){ # if it lies inbetween an exon 
                                    perc <- (start_n-exon_data$V4[exons_orig])/(exon_data$V5[exons_orig]-exon_data$V4[exons_orig])
                                    start_n <- exons$V4[exons_orig]+((exons$V5[exons_orig]-exons$V4[exons_orig])*perc)
                                    inbetween_exon <- "yay"
                                }
                                
                            }
                            
                            if (inbetween_exon == "nay"){
                                for (exons_orig in allrows(exon_data)){ # check first if it is between
                                    if(exons_orig == max(allrows(exon_data))){
                                        
                                    }else{
                                        if(exon_data$V5[exons_orig]<start_n & exon_data$V4[exons_orig+1]>start_n){ # if it lies inbetween an intron 
                                            perc <- (start_n-exon_data$V5[exons_orig])/(exon_data$V4[exons_orig+1]-exon_data$V5[exons_orig])
                                            start_n <- exons$V5[exons_orig]+((exons$V4[exons_orig+1]-exons$V5[exons_orig])*perc)
                                        }
                                    }  
                                } 
                            } 
                            
                            
                        }
                        
                        ###
                        ### look at end value
                        exonic_end <- "nay"
                        inbetween_exon_end <- "nay"
                        for (exons_orig in allrows(exon_data)){ # loop over every exon
                            if(exon_data$V5[exons_orig]==stop_n){ # if one exons has the same start
                                stop_ex <- exons_orig # define this as the start exon
                                exonic_end <- "yay"
                                
                            }
                        } 
                        if (exonic_end == "yay"){
                            stop_n <- exons$V5[stop_ex] # assign the start point from the table of the transformed exons
                        }else{ # the start location is not known :( --> need to find the region to approximate
                            for (exons_orig in allrows(exon_data)){ # check first if it is between 
                                if(exon_data$V4[exons_orig]<stop_n & exon_data$V5[exons_orig]>stop_n){ # if it lies inbetween an exon 
                                    perc <- (stop_n-exon_data$V4[exons_orig])/(exon_data$V5[exons_orig]-exon_data$V4[exons_orig])
                                    stop_n <- exons$V4[exons_orig]+((exons$V5[exons_orig]-exons$V4[exons_orig])*perc)
                                    inbetween_exon_end <- "yay"
                                }
                                
                            }
                            
                            if (inbetween_exon_end == "nay"){
                                for (exons_orig in allrows(exon_data)){ # check first if it is between
                                    if(exons_orig == max(allrows(exon_data))){
                                        
                                    }else{
                                        if(exon_data$V5[exons_orig]<stop_n & exon_data$V4[exons_orig+1]>stop_n){ # if it lies inbetween an intron 
                                            perc <- (stop_n-exon_data$V5[exons_orig])/(exon_data$V4[exons_orig+1]-exon_data$V5[exons_orig])
                                            stop_n <- exons$V5[exons_orig]+((exons$V4[exons_orig+1]-exons$V5[exons_orig])*perc)
                                        }
                                    }
                                } 
                            }
                            
                            
                        }
                        
                        
                        if (position == 1){ # same for first one
                            #max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                            #max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        rect(xleft = max_left+start_n, xright = max_left+stop_n,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                        
                    }else{# if n was a numeric value
                        
                        if (position == 1){ # same for first one
                            max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                        #rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "black",density=15, angle=90, border = FALSE)
                        
                        text(x = max_left+exons$V4[n]+20, y = ref[1]+addtoref[x]-20, labels = n, cex = 0.65)
                        
                        
                        
                    }
                    
                }
            }else{ # if there is no ? in the exon sequence
                for (n in ex_seq){  # loop over exon sequence
                    n <- as.numeric(n)
                    
                    if (n == as.numeric(ex_seq[1])){ # same for first one
                        max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov)#, border = FALSE)
                    #rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "black",density=15, angle=90, border = FALSE)
                    
                    text(x = max_left+exons$V4[n]+20, y = ref[1]+addtoref[x]-20, labels = n, cex = 0.65)
                    
                    
                }
                text(x = max_left - 450*mult, y =  ref[1]+addtoref[x]+25, labels = reference$V1[i]) 
            }
        }
    }
    
    
    )
    
    # add exonic composition of the genetic region 
    
    ex_sorted <- as.data.frame(exons)
    
    addtorefexon <-c(10,27,44,61,78,95,112,129)
    #x <- 7
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                end <- abs(ex_sorted$V5[i])
            }else if((abs(ex_sorted$V4[i])-end)>=15){
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- abs(ex_sorted$V5[i]) 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    
    
    # add legend
    if(mult==as.numeric(-1)){
        xlegend <- 550*mult}else{xlegend <- 7550}
    
    legend(legend = c("< 0.3 FPKM","< 1 FPKM","< 4 FPKM","< 6 FPKM","< 10 FPKM", "> 10 FPKM"), x = xlegend, y = ylim_num+120, fill=c(cols_coverage(6)[1:6]))#, title = "colour legend") 
    
    # add descriptions
    
    if(mult==as.numeric(-1)){
        text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
        
    }else{
        text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        
    }
    
    # add absolute location
    
    chr <- exon_data$V1[1]
    if(mult==as.numeric(-1)){
        gene_start <- min(exon_data$V4)
        gene_end <- max(exon_data$V5)
        gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
        text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
        
        
        
    }else{
        gene_start <- min(exon_data$V4)
        gene_end <- max(exon_data$V5)
        gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
        text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
        text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
        
    }
    
    
    
    # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
    
    text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
    
    
    
}
# two orfs are now displayed correctly in minus stranded genes

PlotAllTranscriptsInExonGraph <- function(exon_data,transcripts,gene_name,vis.ref.file = FALSE, domains = FALSE, references = FALSE){
    
    transcripts <- transcripts[transcripts$V4!="NO ORF FOUND",]
    num_rows <- allrows(transcripts)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    for(i in reps){
        transcripts_reps_rows <- as.numeric(uncollapse(i)[[1]])
        
        PlotTranscriptsInExonGraph(
            exon_data = exons,
            transcripts = transcripts,
            transcript_rows = transcripts_reps_rows,
            gene_name = gene.name,
            vis.ref.file = vis.ref.file,
            domains = domain_df
        )
        
        
    }
    
}
# plots all transcripts in exon graph now only the ones we gave them
GetRefTranscripts <- function(ref_table, exon_df){
    # get reference transcripts from a gtf file and compare it to a selfmade exon table 
    ref <- as.data.frame(ref_table)
    exons <- exon_df
    
    # first remove every row which does not belong to this gene
    # by chromosome & strand first + remove everything which is not transcript & exon
    chr <- exons$V1[1]
    str <- as.character(exons$V7[1])
    ref <- ref[ref$V1==chr,]
    ref <- ref[ref$V7==str,]
    ref <- ref[ref$V3%in%c("exon","transcript"),]
    
    # by location
    min <- min(as.numeric(exons$V4))
    max <- max(as.numeric(exons$V5))
    
    for(i in allrows(ref)){
        if(ref$V4[i]%in%c(min:max)&&ref$V5[i]%in%c(min:max)){
        }else{ref <- ref[-i,]}
    }
    
    # extract transcripts
    tr_table <- c()
    for(i in allrows(ref)){ # for all rows
        
        if(ref$V3[i] == "transcript"){ # if a transcript is described
            tr_name <- c(RemoveElements(uncollapse(ref$V9[i])[[1]][4])) # get tr name
            tr_seq <- c()
            for(n in allrows(ref)){ # check every exon that belong to the transcript and get uniq id
                if(ref$V3[n] == "exon"){
                    #n <- 1
                    ex_tr_name <-  c(RemoveElements(uncollapse(ref$V9[n])[[1]][4]))
                    if(ex_tr_name == tr_name){
                        start <- ref$V4[n]
                        end <- ref$V5[n]
                        exon_checked <- 0
                        for(x in allrows(exons)){
                            if(exons$V4[x]==start&&exons$V5[x]==end){
                                tr_seq <- c(tr_seq,x)     
                            }else{exon_checked <- exon_checked+1}
                        }
                        if(exon_checked == nrow(exons)){
                            tr_seq <- c(tr_seq,"?")
                        }
                    }
                }                    
            }
            
            tr <- c(tr_name, paste(tr_seq, collapse=" "))  
            tr_table <- as.data.frame(rbind(tr_table, tr))
            tr_table <- data.frame(lapply(tr_table, as.character), stringsAsFactors=FALSE)            
        }
    }
    
    
    return(tr_table)   }
# now the input is not a file, but an already loaded file/dataframe
Translate3Frames <- function(seq){
    if(length(seq)==1){
        seq <- s2c(seq)}
    tl1 <- seqinr::translate(seq, frame = 0)
    tl2 <- seqinr::translate(seq, frame = 1)
    tl3 <- seqinr::translate(seq, frame = 2)
    for(n in 1:length(tl1)){
        if(tl1[n]=="*")
        {tl1[n] <- "M"}
    }
    for(n in 1:length(tl2)){
        if(tl2[n]=="*")
        {tl2[n] <- "M"}
    }
    for(n in 1:length(tl3)){
        if(tl3[n]=="*")
        {tl3[n] <- "M"}
    }
    tl1 <- c2s(tl1)
    tl2 <- c2s(tl2)
    tl3 <- c2s(tl3)
    translations <- c(tl1,tl2,tl3)
    return(data.frame(translations, row.names = c("1","2","3"), stringsAsFactors = FALSE))
    
}

PlotExonsWithoutOverlap <- function(vis.ref.file){
    
    
    ex_sorted <- as.data.frame(read.table(file = vis.ref.file, sep = "\t", stringsAsFactors = FALSE))
    
    
    ylim_num <- 1300 # highest point of plot
    
    ref <- c(50+900,100+900) # first line
    addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
    max_left <- 100 # leftmost coordinate
    max_right <- 10100 # rightmost coordinate
    addtorefexon <-c(10,27,44,61,78)
    
    plot(x = c(10200,-300), y = c(0,0), ylim = c(-100,ylim_num), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
    #rect(xleft = 1, xright = 2, ybottom = 0, ytop = 1300, border = "NA", col = "gray")  # plot "y-axis"
    
    rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
    
    for (i in allrows(ex_sorted)){
        
        rect(xleft = max_left+ex_sorted$V4[i], xright = max_left+ex_sorted$V5[i],ybot = ref[1], ytop = ref[2], col = "black",  border = FALSE)
    }
    
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "darkblue", border = FALSE)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                end <- ex_sorted$V5[i]
            }else if((ex_sorted$V4[i]-end)>=15){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = "darkblue", border = FALSE)
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- ex_sorted$V5[i] 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    #text(x = max_left+700, y =  ref[2]+25, labels = "genetic region with transcribed sequence (exons overlapped)", col = "black") 
    
    text(labels = "all sequenced exons without overlap", x = 5100, y = ylim_num -100, cex = 2)
    
}
# added smaller exons on top of genetic region
PlotExonCountDistribution<- function(counts, vis.ref, sample, size.fit = TRUE){
    
    # counts = imported count file from a bedtools result
    # vis.ref from PlotTrInExonGraph
    
    vis.ref <- read.table(file = vis.ref)
    exon_counts <- counts
    
    ### merge to df
    exon_df <- cbind(exon_counts[,c(4,6)],vis.ref[,c(4,5)])
    colnames(exon_df) <- c("exon","count","start","stop")
    
    ### divide count value by exon length (default option)
    if(size.fit == TRUE){
        for (i in allrows(exon_df)){
            exon_df$count[i] <- exon_df$count[i]/(exon_df$stop[i]-exon_df$start[i])
        }
        
    }
    ### fit dataframe
    
    max_val <- max(exon_df$count)
    
    max_val_plot <- plyr::round_any(max_val, accuracy = 100, f=ceiling)
    
    norm_factor <- 1000/max_val_plot
    
    exon_df$count <- exon_df$count*norm_factor
    
    ylim_num <- 1300 # highest point of plot
    
    ref <- c(-200,-100) # first line
    max_left <- 100 # leftmost coordinate
    max_right <- 10100 # rightmost coordinate
    addtorefexon <-c(10,27,44,61,78)
    
    plot(x = c(10200,-300), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "black")
    
    for(i in allrows(exon_df)){
        rect(xleft = max_left+ exon_df$start[i], xright =max_left+ exon_df$stop[i], ybottom = exon_df$count[i]-2, ytop = exon_df$count[i]+2, col = "gray")
    }
    
    
    ex_sorted <- vis.ref
    
    rect(xleft = 0, xright = 0, ybottom = -50, ytop = 1000+100, col = "black")
    
    for (i in allrows(ex_sorted)){
        
        rect(xleft = max_left+ex_sorted$V4[i], xright = max_left+ex_sorted$V5[i],ybot = ref[1], ytop = ref[2], col = "black",  border = FALSE)
    }
    
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                end <- ex_sorted$V5[i]
            }else if((ex_sorted$V4[i]-end)>=15){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- ex_sorted$V5[i] 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    text(labels = max_val_plot, x = -200, y = 1000, cex = 1)
    if(size.fit == TRUE){
    text(labels = "Exon read coverage (per base)", x = 5100, y = ylim_num -100, cex = 2)
    }else{text(labels = "Readcounts of exons", x = 5100, y = ylim_num -100, cex = 2)}
    
    text(labels=as.character(sample),y=1200,x=0, cex = 2, srt = 45)
    
}
# added name and default fitting to exon size (coverage)
TripleLength <- function(char){
    return(length(s2c(char))*3)
}

DoPeptideAndNucleotideSequenceMatch <- function(nucl,pept){
    x <- Translate3Frames(seq = nucl)
    if(x$translations[1]==pept){
        text <- TRUE
    }else{text <- FALSE}
    return(text)
}
# made output Boolean
ReverseComplement <- function(dna_seq){
    rev_comp <-reverseComplement(DNAString(dna_seq))
    return(as.character(rev_comp))
}

ImportReferenceFasta <- function(file, gene_id, translate_to_protein= TRUE){
    fasta <- ReadIn(file)
    
    # line number of different references
    ref_lines <- c()
    for(i in allrows(fasta)){
        if(grepl(">",fasta$V1[i])){
            ref_lines <- c(ref_lines,i)
        }
    }
    
    # num of diff refs
    numbers_ref <- length(ref_lines)
    
    i <- 0
    
    ref_df <- c()
    while(numbers_ref >=1){
        numbers_ref <- numbers_ref-1
        i <- i+1
        
        if(i==length(ref_lines)){
            cols <- c((ref_lines[i]+1):nrow(fasta))
        }else{cols <-c((ref_lines[i]+1):(ref_lines[i+1]-1)) }
        
        # find out whether it is a transcript variant which has published evidence (NR_,NM_) or it is a predicted variant (XM_ or XR_)
        if(grepl(pattern = "NR_", x = as.character(fasta$V1[ref_lines[i]]))){
            evid <- c()
        }else if(grepl(pattern = "NM_", x = as.character(fasta$V1[ref_lines[i]]))){
            evid <- c()
        }else{evid <- "(pred)"}
        
        
        
        ref_df <- rbind(ref_df,c(paste(gene_id,"_ref_",i,evid,collapse='',sep=''),paste(fasta$V1[cols],collapse=""),as.character(fasta$V1[ref_lines[i]])))
        
        
    } 
    
    ref_df <- as.data.frame(ref_df)
    
    ref_df <- taRifx::remove.factors(ref_df)
    
    if(translate_to_protein == TRUE){
        for(n in allrows(ref_df)){
            #i <- 1
            seq_i <- s2c(ref_df$V2[n])
            nucl <- 0       # for safety: count how many possible Nucleotides (ATGC) could be counted
            for (i in seq_i[-length(seq_i)]){
                if(i == "C"){
                    nucl <- nucl+1
                }else if(i == "A"){
                    nucl <- nucl +1
                }else if (i == "T"){
                    nucl <- nucl +1
                }else if (i == "G"){
                    nucl <- nucl+1
                }else if (i == "N"){
                    nucl <- nucl+1
                }
            }
            
            if (nucl==(length(seq_i)-1)){ #if only ATGC + N are in this string --> translate
                ref_df$V2[n] <- c2s(seqinr::translate(s2c(ref_df$V2[n])))
            }
            
            
            
        }
    }
    
    return(ref_df)
    
}
# added translation to prot_seq automatically

AnalyseTranscripts <- function(tr_table, genename_length = 6, relative = FALSE){
    
    # get the affiliation (to a sample) of each transcript
    
    affil <- c()
    type <- c()
    for(i in allrows(tr_table)){
        
        new <- RemoveElements(tr_table$V1[i],start = "first",amount = genename_length+1)
        
        while(s2c(new)[length(s2c(new))]!="."){
            new <- RemoveElements(new,amount=1)
        }
        
        new <- RemoveElements(new,amount=1) # remove the point
        
        affil <- c(affil,new)
        
        while(s2c(new)[length(s2c(new))]!="_"){
            new <- RemoveElements(new,amount=1)
        }
        
        new <- RemoveElements(new,amount=1) # remove the underscore
        
        type <- c(type,new)
        
    }
    
    tr_table <- cbind(tr_table, affil,type) # add it to the transcript
    tr_table <- tr_table[,c(1,9,8,2:7)]
    
    ### change abundance in last coloumn to relative FPKM abundance
    
    affil_num <- length(levels(tr_table$affil)) # how many samples are there?
    
    df_new <- c()
    
    for(i in c(1:affil_num)){ # for each sample
        sample <- levels(tr_table$affil)[i]
        df_tmp <- tr_table[tr_table$affil==sample,]
        
        affil_sum <- 0 # get overall fpkm
        for(n in allrows(df_tmp)){
            affil_sum <- affil_sum + as.numeric(RemoveElements(uncollapse(df_tmp$V7[n])[[1]][4],amount = 1))
        }
        
        for(n in allrows(df_tmp)){ # now calculate ORF-fpkm/overall-FPKM
            rel_fpkm_n <- as.numeric(RemoveElements(uncollapse(df_tmp$V7[n])[[1]][4],amount = 1))/affil_sum
            df_tmp$V7[n] <- paste(df_tmp$V7[n],' rel_FPKM ',round(rel_fpkm_n,digits = 6),';', sep="")
            
        }
        
        df_new <- rbind(df_new,df_tmp)
    }
    
    
    #percentages <-c() # get overall fpkm
    #for(n in allrows(df_new)){
    #    percentages <- c(percentages,as.numeric(RemoveElements(uncollapse(df_new$V7[n])[[1]][8],amount = 1)))
    #}
    
    #plot(percentages)
    
    tr_table <- df_new
    
    ### compose the of the new output dataframe
    
    # absolute or relative values?
    
    if(relative==FALSE){
        fpkm_col <- 4
    }else if(relative==TRUE){
        fpkm_col <- 8
    }
    
    
    # how many types are there?
    type_num <- length(levels(tr_table$type))
    rownames_types <- c()
    for(i in c(1:type_num)){
        # i <- 1
        x <- c(as.character(levels(tr_table$type)[i]),paste(as.character(levels(tr_table$type)[i]),'_meanFPKM',sep='',collapse=''))
        rownames_types <- c(rownames_types,x)
    }
    
    rownames_new <- c("ORF","all","mean_all",c(rownames_types)) 
    # accidentally messed up row and colnames
    # rownames_new should be colnames
    
    rownumber <- length(unique(tr_table$V5))
    
    df_final <- rbind(rownames_new,NA)
    for(i in c(1:(rownumber-1))){
        df_final <- rbind(df_final,NA)
    }
    
    colnames(df_final) <- rownames_new
    df_final <- as.data.frame(df_final[-1,])
    rownames(df_final) <- c(1:length(rownames(df_final)))
    df_final <- remove.factors(df_final)
    
    # find different transcripts 
    x <- length(rownames(tr_table))
    
    tr_table_save <- tr_table
    #tr_table <- tr_table_save
    
    unique_orf_number <- 0
    
    while(x>=1){
        # take the first transcript
        tr_temp <- tr_table[1,]
        unique_orf_number <- unique_orf_number + 1
        # find if there are other orf occurences
        occurences <- c()
        for (i in allrows(tr_table)){
            if(tr_table$V5[i]==tr_temp$V5){
                occurences <- c(occurences,i)
            }
        }
        # make a temporary table
        
        tr_temp <- tr_table[c(occurences),]
        
        # make a summary
        
        ORF <- as.character(tr_temp$V5[1])
        all <- paste(tr_temp$V1,sep=" ",collapse = " ")
        
        mean <- 0
        for(i in allrows(tr_temp)){
            fpkm_i <- RemoveElements(amount = 1,uncollapse(tr_temp$V7[i])[[1]][fpkm_col])
            mean <- mean + as.numeric(fpkm_i)
            
        }
        
        mean <- mean/(length(levels(tr_table_save$affil))) # divide through the number of all samples
        
        df_final$ORF[unique_orf_number] <- ORF
        df_final$all[unique_orf_number] <- all
        df_final$mean_all[unique_orf_number] <- round(as.numeric(mean),digits = 3)
        
        # now a bit harder --> get individual stats for each type (factor in tr_table_save$type)
        for(i in c(1:length(levels(tr_table_save$type)))){
            
            type <- levels(tr_table_save$type)[i]
            
            tr_temp_fac <- tr_temp[tr_temp$type==type,]
            
            if(length(rownames(tr_temp_fac))==0){
                df_final[unique_orf_number,(3+(i*2-1))] <- NA
                df_final[unique_orf_number,(3+(i*2))] <- NA
            }else{
                fac_i <- paste(tr_temp_fac$V1,sep=" ",collapse = " ")
                
                mean_fac <- 0
                
                
                for(n in allrows(tr_temp_fac)){
                    fpkm_n <- RemoveElements(amount = 1,uncollapse(tr_temp_fac$V7[n])[[1]][fpkm_col])
                    mean_fac <- mean_fac + as.numeric(fpkm_n)
                    
                }
                
                # how many samples are in this type
                samples_type <- as.numeric(table(grepl(x = levels(tr_table_save$affil), pattern = type))["TRUE"])
                
                mean_fac <- mean_fac/samples_type # divide through the number of all samples
                
                df_final[unique_orf_number,(3+(i*2-1))] <- fac_i
                df_final[unique_orf_number,(3+(i*2))] <- round(as.numeric(mean_fac),digits = 3)
                
            }# closing from the else statement
        }        
        
        # remove the lines with the same orf
        tr_table <- tr_table[-c(occurences),]
        
        # calculate x - number of lines
        x <- x - length(occurences)
    }
    
    
    
    
    return(df_final)
    
}
# added built-in possibility to give relative or absolute FPKM values

AddSamplesAndRelFPKMToTranscripts <- function(tr_table, genename_length = 6){
    
    # get the affiliation (to a sample) of each transcript
    
    affil <- c()
    type <- c()
    for(i in allrows(tr_table)){
        
        new <- RemoveElements(tr_table$V1[i],start = "first",amount = genename_length+1)
        
        while(s2c(new)[length(s2c(new))]!="."){
            new <- RemoveElements(new,amount=1)
        }
        
        new <- RemoveElements(new,amount=1) # remove the point
        
        affil <- c(affil,new)
        
        while(s2c(new)[length(s2c(new))]!="_"){
            new <- RemoveElements(new,amount=1)
        }
        
        new <- RemoveElements(new,amount=1) # remove the underscore
        
        type <- c(type,new)
        
    }
    
    tr_table <- cbind(tr_table, affil,type) # add it to the transcript
    col_last <- ncol(tr_table)
    tr_table <- tr_table[,c(1,col_last,col_last-1,2:(col_last-2))] # new order
    
    ### change abundance in last coloumn to relative FPKM abundance
    
    affil_num <- length(levels(tr_table$affil)) # how many samples are there?
    
    df_new <- c()
    
    for(i in c(1:affil_num)){ # for each sample
        sample <- levels(tr_table$affil)[i]
        df_tmp <- tr_table[tr_table$affil==sample,]
        
        affil_sum <- 0 # get overall fpkm
        for(n in allrows(df_tmp)){
            affil_sum <- affil_sum + as.numeric(RemoveElements(uncollapse(df_tmp$V7[n])[[1]][4],amount = 1))
        }
        
        for(n in allrows(df_tmp)){ # now calculate ORF-fpkm/overall-FPKM
            rel_fpkm_n <- as.numeric(RemoveElements(uncollapse(df_tmp$V7[n])[[1]][4],amount = 1))/affil_sum
            df_tmp$V7[n] <- paste(df_tmp$V7[n],' rel_FPKM ',round(rel_fpkm_n,digits = 6),';', sep="")
            
        }
        
        df_new <- rbind(df_new,df_tmp)
    }
    
    
    #percentages <-c() # get overall fpkm
    #for(n in allrows(df_new)){
    #    percentages <- c(percentages,as.numeric(RemoveElements(uncollapse(df_new$V7[n])[[1]][8],amount = 1)))
    #}
    
    #plot(percentages)
    
    tr_table <- df_new
    
    return(tr_table)}
# two orfs possible

CompareTranscriptsOverGroups <- function(tr_table, relative = TRUE){
    
    ### first, split the data in order to include the alternative orfs
    tr <- tr_table
    
    tr_table <- c()
    for(i in allrows(tr)){ 
        if(grepl(pattern = "(mo)",x = tr$V1[i])){ # if an alternative ORF exists
            #i <- 8
            line_1 <- tr[i,1:9]
            line_1[1] <- RemoveElements(line_1[1],5)
            line_1[6] <- uncollapse(line_1[6])[[1]][1]
            line_2 <- tr[i,c(1:6,10,8,9)]
            line_2[1] <- RemoveElements(line_2[1],5)
            line_2[1] <- paste(line_2[1],"_alt",sep="",paste="")## add name indication!!!!
            line_2[6] <- uncollapse(line_2[6])[[1]][2]
            
            names(line_2) <-names(line_1)
            rbind(line_1,line_2)
            tr_table <- rbind(tr_table,line_1,line_2)
            
        }else{ # if not simply r-bind the function
            tr_table <- rbind(tr_table,tr[i,-10])
        }        
    }
    
    tr_table <- as.data.frame(tr_table)
    
    if(relative==FALSE){
        fpkm_col <- 4
    }else if(relative==TRUE){
        fpkm_col <- 8
    }
    
    
    # how many types are there?
    type_num <- length(levels(tr_table$type))
    rownames_types <- c()
    for(i in c(1:type_num)){
        # i <- 1
        x <- c(as.character(levels(tr_table$type)[i]),paste(as.character(levels(tr_table$type)[i]),'_meanFPKM',sep='',collapse=''))
        rownames_types <- c(rownames_types,x)
    }
    
    rownames_new <- c("ORF","all","mean_all",c(rownames_types)) 
    # accidentally messed up row and colnames
    # rownames_new should be colnames
    
    rownumber <- length(unique(tr_table$V5))
    
    df_final <- rbind(rownames_new,NA)
    for(i in c(1:(rownumber-1))){
        df_final <- rbind(df_final,NA)
    }
    
    colnames(df_final) <- rownames_new
    df_final <- as.data.frame(df_final[-1,])
    rownames(df_final) <- c(1:length(rownames(df_final)))
    df_final <- remove.factors(df_final)
    
    # find different transcripts 
    x <- length(rownames(tr_table))
    
    tr_table_save <- tr_table
    #tr_table <- tr_table_save
    
    unique_orf_number <- 0
    
    while(x>=1){
        # take the first transcript
        tr_temp <- tr_table[1,]
        unique_orf_number <- unique_orf_number + 1
        # find if there are other orf occurences
        occurences <- c()
        for (i in allrows(tr_table)){
            if(tr_table$V5[i]==tr_temp$V5){
                occurences <- c(occurences,i)
            }
        }
        # make a temporary table
        
        tr_temp <- tr_table[c(occurences),]
        
        # make a summary
        
        ORF <- as.character(tr_temp$V5[1])
        all <- paste(tr_temp$V1,sep=" ",collapse = " ")
        
        mean <- 0
        for(i in allrows(tr_temp)){
            fpkm_i <- RemoveElements(amount = 1,uncollapse(tr_temp$V7[i])[[1]][fpkm_col])
            mean <- mean + as.numeric(fpkm_i)
            
        }
        
        mean <- mean/(length(levels(tr_table_save$affil))) # divide through the number of all samples
        
        df_final$ORF[unique_orf_number] <- ORF
        df_final$all[unique_orf_number] <- all
        df_final$mean_all[unique_orf_number] <- round(as.numeric(mean),digits = 3)
        
        # now a bit harder --> get individual stats for each type (factor in tr_table_save$type)
        for(i in c(1:length(levels(tr_table_save$type)))){
            
            type <- levels(tr_table_save$type)[i]
            
            tr_temp_fac <- tr_temp[tr_temp$type==type,]
            
            if(length(rownames(tr_temp_fac))==0){
                df_final[unique_orf_number,(3+(i*2-1))] <- NA
                df_final[unique_orf_number,(3+(i*2))] <- NA
            }else{
                fac_i <- paste(tr_temp_fac$V1,sep=" ",collapse = " ")
                
                mean_fac <- 0
                
                
                for(n in allrows(tr_temp_fac)){
                    fpkm_n <- RemoveElements(amount = 1,uncollapse(tr_temp_fac$V7[n])[[1]][fpkm_col])
                    mean_fac <- mean_fac + as.numeric(fpkm_n)
                    
                }
                
                # how many samples are in this type
                samples_type <- as.numeric(table(grepl(x = levels(tr_table_save$affil), pattern = type))["TRUE"])
                
                mean_fac <- mean_fac/samples_type # divide through the number of all samples
                
                df_final[unique_orf_number,(3+(i*2-1))] <- fac_i
                df_final[unique_orf_number,(3+(i*2))] <- round(as.numeric(mean_fac),digits = 3)
                
            }# closing from the else statement
        }        
        
        # remove the lines with the same orf
        tr_table <- tr_table[-c(occurences),]
        
        # calculate x - number of lines
        x <- x - length(occurences)
    }
    
    
    
    
    return(df_final)
    
}
 # two orfs possible

CompareReferencesAndORFs <- function(tr_table, ref_table, add_refs = TRUE, add_mass = TRUE, add_length = TRUE){
    
    # tr_table <- output of CompareTranscriptsOverGroups
    # ref_table <- output of ImportReferenceFasta
    
    tr_table <- as.data.frame(cbind(tr_table, NA, NA, NA))
    ref_table <- as.data.frame(ref_table)
    num_col <- length(colnames(tr_table))
    colnames(tr_table)[(num_col-2):(num_col)] <- c("equal","part","include")
    
    # equal
    
    for(i in allrows(ref_table)){
        #i <- 1
        rows_ref <- c(as.numeric(rownames(tr_table[tr_table$ORF==ref_table$V2[i],])))
        for(n in rows_ref){
            if(is.na(tr_table$equal[n])){
                tr_table$equal[n] <- ref_table$V1[i]   
            }else{tr_table$equal[n] <- paste(tr_table$equal[n],ref_table$V1[i],sep=" ",collapse = " ")}
        }
        
    }
    
    # part
    
    for(i in allrows(ref_table)){
        #i <- 1
        # n <- 1
        subject <- ref_table$V2[i]
        for(n in allrows(tr_table)){
            match_res <- matchPattern(pattern = tr_table$ORF[n],subject = subject)
            if(length(start(match_res))>=1){ # if there are matching patterns
                if(is.na(tr_table$part[n])){
                    tr_table$part[n] <- ref_table$V1[i]
                }else{
                    tr_table$part[n] <- paste(tr_table$part[n],ref_table$V1[i],sep=" ", collapse=" ")
                }
                
            }
        }
        
    }
    
    # included
    
    for(i in allrows(ref_table)){
        #i <- 1
        # n <- 1
        pattern <- ref_table$V2[i]
        for(n in allrows(tr_table)){
            match_res <- matchPattern(pattern = pattern, subject = tr_table$ORF[n])
            if(length(start(match_res))>=1){ # if there are matching patterns
                if(is.na(tr_table$include[n])){
                    tr_table$include[n] <- ref_table$V1[i]
                }else{tr_table$include[n] <- paste(tr_table$include[n], ref_table$V1[i],sep = " ", collapse = " " )}
            }
        }
        
    }
    
    #############
    
    if(add_refs == TRUE){ # default
        
        for(i in allrows(ref_table)){
            add_row <- c(ref_table$V2[i],ref_table$V1[i])
            for(i in 3:num_col){
                add_row <- c(add_row, NA)
            }
            # add it to tr_table
            tr_table <- as.data.frame(rbind(tr_table,add_row))
            
        }
        
    }
    
    tr_table <- as.data.frame(cbind(tr_table,NA))
    num_col <- length(colnames(tr_table))
    colnames(tr_table)[num_col] <- "length"
    
    if(add_length == TRUE){
        for(i in allrows(tr_table)){
            tr_table$length[i] <- length(s2c(tr_table$ORF[i]))
        }
    }
    
    
    tr_table <- as.data.frame(cbind(tr_table,NA))
    num_col <- length(colnames(tr_table))
    colnames(tr_table)[num_col] <- "mass_(kDa)"
    
    if(add_length == TRUE){
        for(i in allrows(tr_table)){
            tr_table$`mass_(kDa)`[i] <- CalculatePeptideMass(tr_table$ORF[i],unit = "kDa")
        }
    }
    
    
    return(tr_table)
    
    
}

VerboseInterpretation <- function(compared_transcripts, genename_length = 6){
    
    comp_tr <- compared_transcripts 
    
    # Receptor
    
    rec_name <- uncollapse(comp_tr$all[1])[[1]][1]
    rec_name <- RemoveElements(rec_name,amount = length(s2c(rec_name))-genename_length)
    
    # find out references and non_refernces and seperate
    
    ref_rows <- c()
    for(i in allrows(comp_tr)){
        if(grepl(pattern = "ref", x = comp_tr$all[i])){
            ref_rows <- c(ref_rows,i)
        }
    }
    
    ref_table <- comp_tr[ref_rows,]
    orf_table <- comp_tr[-ref_rows,]
    
    # references
    num_pred <- c()
    for(i in allrows(ref_table)){
        if(grepl(pattern ="pred", x = ref_table$all[i])){
            num_pred <- c(num_pred,i)
        }
    }
    
    num_pred <- length(num_pred)
    num_ref <- length(rownames(ref_table))
    
    # which tissues & samples
    
    num_tissues <- c()
    orf_names <- c()
    for(i in allrows(orf_table)){
        orf_names <- c(orf_names,uncollapse(paste(sep = " ",orf_table$all[i], collapse = ""))[[1]])
    }
    
    overall_num <- length(orf_names)
    
    orf_names_new <- c()
    
    for(i in orf_names){
        orf_names_new <- c(orf_names_new,RemoveElements(char = i,amount = genename_length+1,start = "first"))
    }
    
    sample_names <- c()
    
    for(i in orf_names_new){
        # i <- orf_names_new[1]
        while(s2c(i)[length(s2c(i))]!="."){
            i <- RemoveElements(char = i,amount = 1)
        }
        i <- RemoveElements(char = i,amount = 1)
        
        sample_names<- c(sample_names,i)
        
    }
    
    sample_names <- unique(sample_names)
    
    type_names <- c()
    for(i in sample_names){
        #i = "fat_1"
        while(s2c(i)[length(s2c(i))]!="_"){
            i <- RemoveElements(char = i, amount = 1)
        }
        i <- RemoveElements(char=i,amount=1)
        type_names <- c(type_names,i)
    }
    
    type_names <- unique(type_names)
    
    
    # number of orfs in all types
    
    comp_tr_tmp <- comp_tr[,-c(1:3,length(colnames(comp_tr)),length(colnames(comp_tr))-1,length(colnames(comp_tr))-2,length(colnames(comp_tr))-3,length(colnames(comp_tr))-4)]
    
    reproduced_orfs <- 0
    for(i in allrows(comp_tr_tmp)){
        #i <- 1 
        isitsamplespec <- "no"
        
        for(n in allcols(comp_tr_tmp)){
            if(is.na(comp_tr_tmp[i,n])){
                isitsamplespec <- "yes"
            }
        }
        
        if(isitsamplespec=="no"){
            reproduced_orfs <- c(reproduced_orfs,i)
        }
        
    }
    
    if(length(reproduced_orfs)==1&&reproduced_orfs[1]==0){
        num_reproduced_orfs <- 0
    }else{num_reproduced_orfs <- length(reproduced_orfs)-1}
    
    
    # how many references are there ..> are they found
    
    
    # references found identically
    # references found within
    
    ref_tr_tmp <- comp_tr[,c(length(colnames(comp_tr))-2,length(colnames(comp_tr))-3,length(colnames(comp_tr))-4)]
    
    
    found_refs <- unique(c(ref_tr_tmp$equal,ref_tr_tmp$include))
    
    found_refs <- found_refs[!is.na(found_refs)]    
    
    f_refs <- c()
    for(i in found_refs){
        f_refs <- c(f_refs,uncollapse(i)[[1]])
    }
    
    f_refs <- unique(f_refs)
    
    # verbose output
    print(paste("Gene: ",rec_name,'. Tissues: ',paste(type_names,sep=" & ", collapse=" & "),'. Samples: ',paste(sample_names,sep=" & ", collapse=" & "),".",sep=""))
    print(paste("Overall Transcripts: ",overall_num,'. Different/unique transcripts: ',length(rownames(orf_table)),'. Number of ORFs found in all tissues: ',num_reproduced_orfs,'.',sep=""))
    print(paste("Amount of annotated references: ",num_ref," (of which ",num_pred," are predicted.) ",length(f_refs),' references were found in the transcripts: ',paste(f_refs,sep=" & ",collapse=" & "),'.',sep=""))
}
# small mistakes were bettered

ComparativeTranscriptPlot_alt <- function(exon_data, excluded_exons = FALSE, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(allrows(comp_tr))
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    ### extract information about each compared transcript from transcript_list and add it to comp_tr
    
    exon_number_seq <- c()
    dna_seq <- c()
    orf_seq <- c()
    
    for(i in allrows(comp_tr)){
        # i <- 6        
        example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
        
        if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
            
            example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
        }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }else{
            example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
        }
        
    }
    
    # add to comp_tr 
    
    comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
    
    ### exon preprocessing
    
    # remove excluded_exons from ex_exons
    
    suppressWarnings(if(excluded_exons==FALSE){
        exons <- exon_df
        excl_exons <- c(0)
        excluded_exons <- c(0)
    }else{ # remove the exons but before: check if they are needed in either a domain or a transcript
        needed_exons <- c()
        if(exist.domains=="yes"){
            for(n in allrows(domains)){
                needed_exons <- c(needed_exons,uncollapse(domains$domain_locations[n])[[1]])
                needed_exons <- unique(floor(as.numeric(needed_exons)))
            }
        }
        
        for(n in allrows(comp_tr)){
            needed_exons <- unique(as.numeric(c(needed_exons,uncollapse(comp_tr$exon_number_seq[n])[[1]])))
        }
        
        excl_exons <- c()
        for(n in excluded_exons){
            if(n%in%needed_exons){
                print("Some exons could not be omitted from location fitting.")
            }else{excl_exons <- as.numeric(c(excl_exons,n))}
        }
        
        exons <- exon_df[-excl_exons,]
    }
    )
    
    # location processing 
    # find first exon beginning and substract this value
    
    min_start <- min(exons$V4)
    exons$V4 <- exons$V4 - min_start
    exons$V5 <- exons$V5 - min_start
    
    # find overlap groups
    
    # for this, the rownames are important, but the rownames are also important for later re-adding the excluded exons
    
    rownames_save <- rownames(exons)
    
    rownames(exons) <- allrows(exons)
    
    FindFirstOverlapGroup <- function(exon_df, rows){
        df_ex <- as.data.frame(exon_df)
        for (i in c(rows)){   # for every exon
            if (i == min(rows)){    # if it is the first exon of the input
                group <- c(i)
                range_a <- c(df_ex$V4[i]:df_ex$V5[i])   # define first range
            }else{      # for every other exon
                range_tmp <- c(df_ex$V4[i]:df_ex$V5[i]) # define its own range (_tmp)    
                is_part_of_a <- 0         # set this parameter to 0
                for (n in range_tmp){   # for ever value in range_tmp
                    if(n%in%range_a){   # check if it is in range_a
                        is_part_of_a <- is_part_of_a+1  # and if this is the case, increment is_part_of by one 
                    }
                }
                if (is_part_of_a >= 1){   # if is_part_of >= 1, this means, that range_temp is part of range_a
                    range_a <- unique(c(range_a, range_tmp)) # merge them 
                    group <- c(group,i)   # and add to group
                } else{
                    # not_in_group <- c(not_in_group,i)
                }
            }
        }
        #print(group)
        return(group)
    }
    
    # group <- 18:21
    firsts <- c(0) # make sure to delete the zero later --> zero gets removed somehow automatically, yeah!
    rowsmax <- length(rownames(exons))
    longest <- c(0)
    #i <- 2
    for (i in c(1:rowsmax)){ # for every row, (theoretically, it ends before, but it COULD be the case, that there are rowmax different exons with no overlap overall)
        if (max(firsts)<length(rownames(exons))){#if the highest/last exon added to first is still lower than the last exon --> do the analysis
            if (i==1){
                group <- FindFirstOverlapGroup(exons,rows = (c(1:rowsmax)))
                firsts <- group[length(FindFirstOverlapGroup(exons,rows = (c(1:rowsmax))))] # get the last of the first overlap group
                group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                longest <- rownames(group_tmp[which.max(group_tmp$V5),])[1] # if there are multiple --> take the first, it doesnt matter                
                
            }else{
                group <- FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax)))
                firsts <- c(firsts,group[length(FindFirstOverlapGroup(exons,rows = (c(max(firsts+1):rowsmax))))]) # get the last exon of the next overlap group
                group_tmp <- exons[group,] # make a group and find out the exon with the highest end 
                longest <- as.numeric(c(longest,rownames(group_tmp[which.max(group_tmp$V5),])[1]))                
                
            }
        }}
    # its okay if there is an error message
    
    #firsts # these are all the last exons of the overlap group
    firsts <- firsts[-length(firsts)] # remove the last one
    seconds <- firsts+1 # these are all exons whose are the first exons of the following overlap group
    longest <- as.numeric(longest[-length(longest)])
    
    exon_dist_rows <- as.data.frame(cbind(firsts,seconds,longest)) # merge to dataframe
    
    subtrExonDist <- function(exon_table,before,after,longest){ # make function
        subtr <- exon_table$V4[after]-exon_table$V5[longest] - 100   # get intron length - 100
        exon_table$V4[c(after:length(rownames(exon_table)))] <- exon_table$V4[c(after:length(rownames(exon_table)))]-subtr # subtract this value
        exon_table$V5[c(after:length(rownames(exon_table)))] <- exon_table$V5[c(after:length(rownames(exon_table)))]-subtr # subtract this value 
        return(exon_table) # return the new exon table
    }
    
    #i <- 1
    for (i in c(1:length(rownames(exon_dist_rows)))){ # for every intron beginning
        bf <- exon_dist_rows$firsts[i]  # get start number
        af <- exon_dist_rows$seconds[i] # get end number
        lon <- exon_dist_rows$longest[i]
        af_minus_long <- as.numeric(exons$V4[af])-as.numeric(exons$V5[lon]) # check distance/intron length
        if (as.numeric(af_minus_long) >= 100){ # if bigger or equal 100
            exons <- subtrExonDist(exon_table = exons, before = bf, after = af, longest = lon) # transform distance to 100
        }
        
    }
    
    
    # normalize for possible place: 10000 Units
    
    max_end <- max(exons$V5) # get highest
    
    fact <- 10000/max_end # get multiplicator
    exons$V4 <- fact*exons$V4 # multiplicate
    exons$V5 <- fact*exons$V5
    exons$V4 <- round(digits = 1,x = exons$V4) # round the results
    exons$V5 <- round(digits = 1, x = exons$V5)
    
    exons$V4 <- exons$V4*mult
    exons$V5 <- exons$V5*mult
    
    rownames(exons) <- rownames_save # assign rownames back
    
    exons_fitted_save <- exons # save the fittet exons
    
    # re-add the excluded exons to the dataframe (in order to correctly plot the other exons)
    
    excl_exons_save <- excl_exons 
    excl_exons <- sort(excl_exons)
    if(excluded_exons[1]==0){
        excl_exons_save <- c(0)
    }else{
        
        if(excl_exons[1]==1){
            merged_exons <- exon_df[excl_exons[1],]
            excl_exons <- excl_exons[-1]
        }else{merged_exons <- c()}
        
        x <- 1
        while(length(excl_exons)>=1){
            if(x>length(rownames(exons))){
                merged_exons <- rbind(merged_exons,exon_df[c(excl_exons),])  
                excl_exons <- c()
            }else{
                if(as.numeric(rownames(exons[x,]))<excl_exons[1]){
                    merged_exons <- rbind(merged_exons,exons[x,])
                    x <- x + 1
                }else{
                    merged_exons <- rbind(merged_exons,exon_df[excl_exons[1],])
                    excl_exons <- excl_exons[-1]
                }
            } 
            
        }        
        
        exons <- merged_exons
        
    }
    
    
    ### plot: basic parameters (don't forget minus strand)
    
    ylim_num <- 1300 # highest point of plot
    
    ref <- c(50+900,100+900) # first line
    addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
    max_left <- 100 * mult # leftmost coordinate
    max_right <- 10100 * mult# rightmost coordinate
    col_ref <- "darkblue"    
    
    if(mult==as.numeric(-1)){
        additional_place <- 200 # important for heatmap - develop that for - strand genes later!
        additional_place_2 <- 1200
    }else{additional_place <- 1200
    additional_place_2 <- 0}
    
    plot_x1 <- (additional_place+10200) * mult
    plot_x2 <- (-300-additional_place_2) * mult
    plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
    
    
    
    ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
    
    rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
    
    # plot every possible exon in the highest line
    
    for (i in allrows(exons_fitted_save)){
        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
    }
    
    #   add exonic composition of the genetic region 
    
    ex_sorted <- as.data.frame(exons_fitted_save)
    
    addtorefexon <-c(10,27,44,61,78,95,112,129)
    
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                end <- abs(ex_sorted$V5[i])
            }else if((abs(ex_sorted$V4[i])-end)>=15){
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- abs(ex_sorted$V5[i]) 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    # add descriptions
    
    if(mult==as.numeric(-1)){
        text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
        
    }else{
        text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        
    }
    
    # plot: domains
    suppressWarnings(if(domains == FALSE){
        # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
    }else{
        
        dom <- domains
        
        for(i in allrows(dom)){
            ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
            
            
            firstexon <- floor(ex_seq[1]) # get integer exon
            perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
            left_coord <- exons$V4[firstexon]+((exons$V5[firstexon]-exons$V4[firstexon])*perc_f_exon)
            
            
            lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
            perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
            if(ex_seq[length(ex_seq)]%%1 > 0){
                right_coord <- exons$V4[lastexon]+((exons$V5[lastexon]-exons$V4[lastexon])*perc_l_exon)
            }else{right_coord <- exons$V5[lastexon]}
            
            ex_seq[1] <- firstexon
            ex_seq[length(ex_seq)] <- lastexon
            if(length(unique(ex_seq))==1){
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }else{for(n in ex_seq){
                if(n==firstexon){
                    rect(xleft = max_left+left_coord, xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else if(n==lastexon){
                    rect(xleft = max_left+exons$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else{rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }
                
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }}
            
            
            rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
            
        }
    })
    
    
    
    
    # add absolute location
    
    chr <- exon_data$V1[1]
    suppressWarnings(if(excluded_exons == 0){
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
    }else{
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[-c(excl_exons_save),4])
            gene_end <- max(exon_data[-c(excl_exons_save),5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[-c(excl_exons_save),4])
            gene_end <- max(exon_data[-c(excl_exons_save),5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
    }
    )
    
    
    
    
    # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
    
    text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
    
    ### plot: comp_tr
    
    # basic color
    
    col_cov_utr <- "orangered2"
    col_cov_orf <- "darkred"
    
    # careful with first & last exon
    # n <- 1
    
    for (n in allrows(comp_tr)){  # for chosen transcript rows
        ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        # make x == n linenumber of this transcript
        
        x <- n
        
        # get exons which are in the ORF and in the UTR
        
        sequence <- as.character(comp_tr$dna_seq[n])
        orf <- as.character(comp_tr$orf_seq[n])
        
        orf_position <- matchPattern(orf, sequence) # find orf in sequence
        
        orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
        orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
        
        exon_sum <- 0
        for (i in ex_seq){ # get exon sum 
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        }
        
        orf_start <- orf_start_rel * exon_sum
        orf_end <- orf_end_rel * exon_sum
        
        exon_sum_max <- exon_sum
        exon_sum <- 0
        orf_exon_first <- NA
        orf_exon_last <- NA
        UTR_exons <- c()
        
        for(i in ex_seq){ # find exon in which the orf begins
            #i <- 190
            i <- as.numeric(i)
            if(is.na(orf_exon_first)){
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum<=orf_start){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_first)){
                    orf_exon_first <- i
                }}
        }
        
        ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
        
        exon_sum <- exon_sum_max
        for(i in ex_seq_back){
            i <- as.numeric(i)
            exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
            if(exon_sum>=orf_end){
                UTR_exons <- c(UTR_exons,i)
            }else if(is.na(orf_exon_last)){
                orf_exon_last <- i
            }
        }
        
        
        # group all exons which are not (or not completely) part of the orf  
        
        UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
        
        # special case: the two exons in which the orf starts/end - find position within exon
        
        # orf start
        
        ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
        
        exon_sum <- 0
        for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        }
        
        sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
        orf_start_exon <- orf_start-sum_before # subtract this from the orf start
        # now we have the value within the first exon
        
        
        if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
            orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
            orf_start_exon_end <- exons$V4[orf_exon_first]      
        }else{
            orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
            orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
        }
        
        # now for orf end
        
        ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
        
        exon_sum <- 0
        for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        } 
        
        sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
        orf_end_exon <- orf_end-sum_before # subtract this from the orf end
        # now we have the value within the last exon
        
        if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
            orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
            orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
        }else{
            orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
            orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
        }
        
        
        # plot exons
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        if(ORFS_only==TRUE){
            
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                    
                    max_left_ex <- 10+max_left+exons$V4[orf_exon_first]
                    max_right_ex <- max_left+exons$V5[orf_exon_last]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    
                    
                }
                
                if(i %in% c(orf_exon_last,orf_exon_first)){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                    text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    
                }else if(i%in%UTR_exons){
                    #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                    
                }else{
                    
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                }
                
                
                
                if(orf_exon_first==orf_exon_last){
                    rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }}
                
                
            }
        }else{
            
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){
                    rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }}
                
                
            }
        } #closing from if ORFS only else 
        #if(mult==as.numeric(-1)){
        #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
        #}else{
        #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
        #}
    }
    
    ### add IDs for each transcript
    
    suppressWarnings(if(IDs == FALSE){
        # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
        IDs <- c("a","b","","d","e","f","g","")
    }else{
        if(mult==as.numeric(-1)){
            for(i in c(1:length(IDs))){
                
                text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
            }
        }else{
            for(i in c(1:length(IDs))){
                text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
            }
            
        }
        
    })
    
    
    
    
    ### plot: heatmap 
    
    if(mult==as.numeric(1)){
        # define the room
        
        #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
        
        # divide it by the number of differents tissues
        
        constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
        tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
        
        tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
        num_tissues <- length(tissue_cols)
        tissues <- colnames(tissue_df)[tissue_cols]
        tissue_place <- additional_place/num_tissues
        for(i in c(1:num_tissues)){
            text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
        }
        
        if(FPKM == TRUE){
            text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
            
            col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    fpkm_val <- tissue_df[i,(2*n)]
                    if(is.na(fpkm_val)){
                        fpkm_col <- "white"
                    }else if(fpkm_val<=0.1){
                        fpkm_col <- col_table_fpkm(7)[1]
                    }else if(fpkm_val<=0.5){
                        fpkm_col <- col_table_fpkm(7)[2]
                    }else if(fpkm_val<=0.999){
                        fpkm_col <- col_table_fpkm(7)[3]
                    }else if(fpkm_val<=1.5){
                        fpkm_col <- col_table_fpkm(7)[4]
                    }else if(fpkm_val<=2.5){
                        fpkm_col <- col_table_fpkm(7)[5]
                    }else if(fpkm_val<=5){
                        fpkm_col <- col_table_fpkm(7)[6]
                    }else if(fpkm_val>5){
                        fpkm_col <- col_table_fpkm(7)[7]
                    }
                    
                    rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                }
            }
            # legend below
            legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
            
        }else{
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                    }
                    
                    rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
        }
    }else if(mult==as.numeric(-1)){ # if mult == -1
        
        # define the room
        
        # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
        
        # divide it by the number of differents tissues
        
        constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
        tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
        
        tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
        num_tissues <- length(tissue_cols)
        tissues <- colnames(tissue_df)[tissue_cols]
        tissue_place <- additional_place_2/num_tissues
        for(i in c(1:num_tissues)){
            text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
        }
        
        if(FPKM == TRUE){
            text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
            
            col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    fpkm_val <- tissue_df[i,(2*n)]
                    if(is.na(fpkm_val)){
                        fpkm_col <- "white"
                    }else if(fpkm_val<=0.1){
                        fpkm_col <- col_table_fpkm(7)[1]
                    }else if(fpkm_val<=0.5){
                        fpkm_col <- col_table_fpkm(7)[2]
                    }else if(fpkm_val<=0.999){
                        fpkm_col <- col_table_fpkm(7)[3]
                    }else if(fpkm_val<=1.5){
                        fpkm_col <- col_table_fpkm(7)[4]
                    }else if(fpkm_val<=2.5){
                        fpkm_col <- col_table_fpkm(7)[5]
                    }else if(fpkm_val<=5){
                        fpkm_col <- col_table_fpkm(7)[6]
                    }else if(fpkm_val>5){
                        fpkm_col <- col_table_fpkm(7)[7]
                    }
                    
                    rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                }
            }
            # legend below
            legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
            
        }else{
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                    }
                    
                    rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
        }
        
        
    }
    
    ### legend
    if(mult==as.numeric(-1)){
        legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
    }else{
        legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        
    }
    
}
# removed non-coding exons

ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    ### extract information about each compared transcript from transcript_list and add it to comp_tr
    
    exon_number_seq <- c()
    dna_seq <- c()
    orf_seq <- c()
    
    for(i in allrows(comp_tr)){
        # i <- 6        
        example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
        
        if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
            
            example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
        }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }else{
            example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
        }
        
    }
    
    # add to comp_tr 
    
    comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
    
    ### exon preprocessing
    
    # load vis ref file
    exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
    exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
    ### plot: basic parameters (don't forget minus strand)
    
    ylim_num <- 1300 # highest point of plot
    ref <- c(50+900,100+900) # first line
    addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
    max_left <- 100 * mult # leftmost coordinate
    max_right <- 10100 * mult# rightmost coordinate
    col_ref <- "darkblue"    
    
    if(mult==as.numeric(-1)){
        additional_place <- 200 # important for heatmap - develop that for - strand genes later!
        additional_place_2 <- 1200
    }else{additional_place <- 1200
    additional_place_2 <- 0}
    
    plot_x1 <- (additional_place+10200) * mult
    plot_x2 <- (-300-additional_place_2) * mult
    plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
    
    
    
    ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
    
    rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
    
    # plot every possible exon in the highest line
    
    for (i in allrows(exons_fitted_save)){
        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
    }
    
    #   add exonic composition of the genetic region 
    
    ex_sorted <- as.data.frame(exons_fitted_save)
    
    addtorefexon <-c(10,27,44,61,78,95,112,129)
    
    keep <- allrows(ex_sorted)
    x <- 1
    while(length(keep)>=1){
        new_keep <- c()
        for(i in keep){
            if(i == keep[1]){
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                end <- abs(ex_sorted$V5[i])
            }else if((abs(ex_sorted$V4[i])-end)>=15){
                #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                end <- abs(ex_sorted$V5[i]) 
            }else{new_keep <- c(new_keep,i)}
            
        }
        end <- c()
        keep <- new_keep
        x <- x +1
    }
    
    # add descriptions
    
    if(mult==as.numeric(-1)){
        text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
        
    }else{
        text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
        text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
        
    }
    
    # plot: domains
    suppressWarnings(if(domains == FALSE){
        # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
    }else{
        
        dom <- domains
        
        for(i in allrows(dom)){
            ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
            
            
            firstexon <- floor(ex_seq[1]) # get integer exon
            perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
            left_coord <- exons$V4[firstexon]+((exons$V5[firstexon]-exons$V4[firstexon])*perc_f_exon)
            
            
            lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
            perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
            if(ex_seq[length(ex_seq)]%%1 > 0){
                right_coord <- exons$V4[lastexon]+((exons$V5[lastexon]-exons$V4[lastexon])*perc_l_exon)
            }else{right_coord <- exons$V5[lastexon]}
            
            ex_seq[1] <- firstexon
            ex_seq[length(ex_seq)] <- lastexon
            if(length(unique(ex_seq))==1){
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }else{for(n in ex_seq){
                if(n==firstexon){
                    rect(xleft = max_left+left_coord, xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else if(n==lastexon){
                    rect(xleft = max_left+exons$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }else{rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                }
                
                mean_dom <- max_left + (left_coord+right_coord)/2
                text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                
            }}
            
            
            rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
            
        }
    })
    
    
    
    
    # add absolute location
    
    chr <- exon_data$V1[1]
    suppressWarnings(if(excluded_exons == 0){
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
    }else{
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[-c(excl_exons_save),4])
            gene_end <- max(exon_data[-c(excl_exons_save),5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[-c(excl_exons_save),4])
            gene_end <- max(exon_data[-c(excl_exons_save),5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
    }
    )
    
    
    
    
    # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
    
    text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
    
    ### plot: comp_tr
    
    # basic color
    
    col_cov_utr <- "orangered2"
    col_cov_orf <- "darkred"
    
    # careful with first & last exon
    # n <- 1
    
    for (n in allrows(comp_tr)){  # for chosen transcript rows
        ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        # make x == n linenumber of this transcript
        
        x <- n
        
        # get exons which are in the ORF and in the UTR
        
        sequence <- as.character(comp_tr$dna_seq[n])
        orf <- as.character(comp_tr$orf_seq[n])
        
        orf_position <- matchPattern(orf, sequence) # find orf in sequence
        
        orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
        orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
        
        exon_sum <- 0
        for (i in ex_seq){ # get exon sum 
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        }
        
        orf_start <- orf_start_rel * exon_sum
        orf_end <- orf_end_rel * exon_sum
        
        exon_sum_max <- exon_sum
        exon_sum <- 0
        orf_exon_first <- NA
        orf_exon_last <- NA
        UTR_exons <- c()
        
        for(i in ex_seq){ # find exon in which the orf begins
            #i <- 190
            i <- as.numeric(i)
            if(is.na(orf_exon_first)){
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum<=orf_start){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_first)){
                    orf_exon_first <- i
                }}
        }
        
        ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
        
        exon_sum <- exon_sum_max
        for(i in ex_seq_back){
            i <- as.numeric(i)
            exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
            if(exon_sum>=orf_end){
                UTR_exons <- c(UTR_exons,i)
            }else if(is.na(orf_exon_last)){
                orf_exon_last <- i
            }
        }
        
        
        # group all exons which are not (or not completely) part of the orf  
        
        UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
        
        # special case: the two exons in which the orf starts/end - find position within exon
        
        # orf start
        
        ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
        
        exon_sum <- 0
        for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        }
        
        sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
        orf_start_exon <- orf_start-sum_before # subtract this from the orf start
        # now we have the value within the first exon
        
        
        if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
            orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
            orf_start_exon_end <- exons$V4[orf_exon_first]      
        }else{
            orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
            orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
        }
        
        # now for orf end
        
        ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
        
        exon_sum <- 0
        for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
            i <- as.numeric(i)
            exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
        } 
        
        sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
        orf_end_exon <- orf_end-sum_before # subtract this from the orf end
        # now we have the value within the last exon
        
        if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
            orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
            orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
        }else{
            orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
            orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
        }
        
        
        # plot exons
        if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
            ex_seq <- ex_seq[length(ex_seq):1]
        }
        
        if(ORFS_only==TRUE){
            
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                    
                    max_left_ex <- 10+max_left+exons$V4[orf_exon_first]
                    max_right_ex <- max_left+exons$V5[orf_exon_last]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    
                    
                }
                
                if(i %in% c(orf_exon_last,orf_exon_first)){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                    text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    
                }else if(i%in%UTR_exons){
                    #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                    
                }else{
                    
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                }
                
                
                
                if(orf_exon_first==orf_exon_last){
                    rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }}
                
                
            }
        }else{
            
            for (i in ex_seq){  # loop over exon sequence
                i <- as.numeric(i)
                
                if (i == min(as.numeric(ex_seq))){ # plot exon line
                    max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                    max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                    rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                    #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                    
                }
                
                if(i%in%UTR_exons){
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                }else{
                    rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                }
                
                text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                
                if(orf_exon_first==orf_exon_last){
                    rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    
                }else{
                    if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                        rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }
                    
                    if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                        rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }}
                
                
            }
        } #closing from if ORFS only else 
        #if(mult==as.numeric(-1)){
        #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
        #}else{
        #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
        #}
    }
    
    ### add IDs for each transcript
    
    suppressWarnings(if(IDs == FALSE){
        # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
        
        #IDs <- c("a","b","","d","e","f","g","")
    }else{
        if(mult==as.numeric(-1)){
            for(i in c(1:length(IDs))){
                
                text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
            }
        }else{
            for(i in c(1:length(IDs))){
                text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
            }
            
        }
        
    })
    
    
    
    
    ### plot: heatmap 
    
    if(mult==as.numeric(1)){
        # define the room
        
        #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
        
        # divide it by the number of differents tissues
        
        constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
        tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
        
        tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
        num_tissues <- length(tissue_cols)
        tissues <- colnames(tissue_df)[tissue_cols]
        tissue_place <- additional_place/num_tissues
        for(i in c(1:num_tissues)){
            text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
        }
        
        if(FPKM == TRUE){
            text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
            
            col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    fpkm_val <- tissue_df[i,(2*n)]
                    if(is.na(fpkm_val)){
                        fpkm_col <- "white"
                    }else if(fpkm_val<=0.1){
                        fpkm_col <- col_table_fpkm(7)[1]
                    }else if(fpkm_val<=0.5){
                        fpkm_col <- col_table_fpkm(7)[2]
                    }else if(fpkm_val<=0.999){
                        fpkm_col <- col_table_fpkm(7)[3]
                    }else if(fpkm_val<=1.5){
                        fpkm_col <- col_table_fpkm(7)[4]
                    }else if(fpkm_val<=2.5){
                        fpkm_col <- col_table_fpkm(7)[5]
                    }else if(fpkm_val<=5){
                        fpkm_col <- col_table_fpkm(7)[6]
                    }else if(fpkm_val>5){
                        fpkm_col <- col_table_fpkm(7)[7]
                    }
                    
                    rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                }
            }
            # legend below
            legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
            
        }else{
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                    }
                    
                    rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
        }
    }else if(mult==as.numeric(-1)){ # if mult == -1
        
        # define the room
        
        # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
        
        # divide it by the number of differents tissues
        
        constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
        tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
        
        tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
        num_tissues <- length(tissue_cols)
        tissues <- colnames(tissue_df)[tissue_cols]
        tissue_place <- additional_place_2/num_tissues
        for(i in c(1:num_tissues)){
            text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
        }
        
        if(FPKM == TRUE){
            text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
            
            col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    fpkm_val <- tissue_df[i,(2*n)]
                    if(is.na(fpkm_val)){
                        fpkm_col <- "white"
                    }else if(fpkm_val<=0.1){
                        fpkm_col <- col_table_fpkm(7)[1]
                    }else if(fpkm_val<=0.5){
                        fpkm_col <- col_table_fpkm(7)[2]
                    }else if(fpkm_val<=0.999){
                        fpkm_col <- col_table_fpkm(7)[3]
                    }else if(fpkm_val<=1.5){
                        fpkm_col <- col_table_fpkm(7)[4]
                    }else if(fpkm_val<=2.5){
                        fpkm_col <- col_table_fpkm(7)[5]
                    }else if(fpkm_val<=5){
                        fpkm_col <- col_table_fpkm(7)[6]
                    }else if(fpkm_val>5){
                        fpkm_col <- col_table_fpkm(7)[7]
                    }
                    
                    rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                }
            }
            # legend below
            legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
            
        }else{
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                    }
                    
                    rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
        }
        
        
    }
    
    ### legend
    if(mult==as.numeric(-1)){
        legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
    }else{
        legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        
    }
    
}
# added vis.ref.file as a must & changed strange behaviour regarding orf coloring
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = tr_compared$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[6] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- compared_transcripts[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 6        
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons$V4[firstexon]+((exons$V5[firstexon]-exons$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons$V4[lastexon]+((exons$V5[lastexon]-exons$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons$V4[n], xright = max_left+exons$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "orangered2"
        col_cov_orf <- "darkred"
        
        # careful with first & last exon
        # n <- 1
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # get exons which are in the ORF and in the UTR
            
            sequence <- as.character(comp_tr$dna_seq[n])
            orf <- as.character(comp_tr$orf_seq[n])
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            if(ORFS_only==TRUE){
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        
                        max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                        
                        max_left_ex <- 10+max_left+exons$V4[orf_exon_first]
                        max_right_ex <- max_left+exons$V5[orf_exon_last]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        
                        
                    }
                    
                    if(i %in% c(orf_exon_last,orf_exon_first)){
                        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                        text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                    }else if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                        
                    }else{
                        
                        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
            }else{
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
            } #closing from if ORFS only else 
            #if(mult==as.numeric(-1)){
            #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            #}else{
            #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            #}
        }
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# now: plots every transcripts (without references indicated by "ref") which is given in compared_transcripts
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = tr_compared$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[6] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- compared_transcripts[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 6        
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "orangered2"
        col_cov_orf <- "darkred"
        
        # careful with first & last exon
        # n <- 1
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # get exons which are in the ORF and in the UTR
            
            sequence <- as.character(comp_tr$dna_seq[n])
            orf <- as.character(comp_tr$orf_seq[n])
            
            orf_position <- matchPattern(orf, sequence) # find orf in sequence
            
            orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
            orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
            
            exon_sum <- 0
            for (i in ex_seq){ # get exon sum 
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            orf_start <- orf_start_rel * exon_sum
            orf_end <- orf_end_rel * exon_sum
            
            exon_sum_max <- exon_sum
            exon_sum <- 0
            orf_exon_first <- NA
            orf_exon_last <- NA
            UTR_exons <- c()
            
            for(i in ex_seq){ # find exon in which the orf begins
                #i <- 190
                i <- as.numeric(i)
                if(is.na(orf_exon_first)){
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum<=orf_start){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_first)){
                        orf_exon_first <- i
                    }}
            }
            
            ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
            
            exon_sum <- exon_sum_max
            for(i in ex_seq_back){
                i <- as.numeric(i)
                exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                if(exon_sum>=orf_end){
                    UTR_exons <- c(UTR_exons,i)
                }else if(is.na(orf_exon_last)){
                    orf_exon_last <- i
                }
            }
            
            
            # group all exons which are not (or not completely) part of the orf  
            
            UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
            
            # special case: the two exons in which the orf starts/end - find position within exon
            
            # orf start
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            }
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
            orf_start_exon <- orf_start-sum_before # subtract this from the orf start
            # now we have the value within the first exon
            
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                orf_start_exon_end <- exons$V4[orf_exon_first]      
            }else{
                orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
            }
            
            # now for orf end
            
            ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
            
            exon_sum <- 0
            for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                i <- as.numeric(i)
                exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
            } 
            
            sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
            orf_end_exon <- orf_end-sum_before # subtract this from the orf end
            # now we have the value within the last exon
            
            if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
            }else{
                orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
            }
            
            
            # plot exons
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            if(ORFS_only==TRUE){
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                        
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                        max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        
                        
                    }
                    
                    if(i %in% c(orf_exon_last,orf_exon_first)){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                    }else if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                        
                    }else{
                        
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
            }else{
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
            } #closing from if ORFS only else 
            #if(mult==as.numeric(-1)){
            #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
            #}else{
            #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
            #}
        }
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# adapted it to minus stranded genes
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = tr_compared$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- compared_transcripts[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "orangered2"
        col_cov_orf <- "darkred"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# fixed error where NA was being in the comp_tr
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = tr_compared$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- compared_transcripts[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        col_table_perc <- colorRampPalette(colors = c('darkred',"grey",'darkblue'))
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('darkblue','white','darkred'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                            rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        }
                        
                        
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                            rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                            
                        }
                        
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# changed colors in heatmap and transcripts
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = tr_compared$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- compared_transcripts[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        col_table_perc <- colorRampPalette(colors = c('grey',"darkred",'darkblue'))
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('darkblue','white','darkred'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                            rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        }
                        
                        
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                            rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                            
                        }
                        
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# changed colors again
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = compared_transcripts$all)!=TRUE,]
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        col_table_perc <- colorRampPalette(colors = c('grey',"darkred",'darkblue'))
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('darkblue','white','darkred'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                            rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        }
                        
                        
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                            rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                            
                        }
                        
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# small adjustments, now no ORF is being left behind
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = compared_transcripts$all)!=TRUE,]
    comp_tr <- comp_tr[comp_tr$ORF!="NO ORF FOUND",]
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        col_table_perc <- colorRampPalette(colors = c('grey',"darkred",'darkblue'))
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('darkblue','white','darkred'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                            rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        }
                        
                        
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                            rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                            
                        }
                        
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# included comp_tr_save
ComparativeTranscriptPlot <- function(exon_data, transcript_list, compared_transcripts, FPKM = FALSE, vis.ref.file = FALSE, domains = FALSE, IDs = FALSE, ORFS_only = FALSE){
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts[grepl(pattern = "ref",x = compared_transcripts$all)!=TRUE,]
    comp_tr <- comp_tr[comp_tr$ORF!="NO ORF FOUND",]
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            example_tr <- uncollapse(comp_tr$all[i])[[1]][1]
            
            if(grepl("alt",x = example_tr)){ # if the alternative orf is needed
                
                example_tr <- paste(RemoveElements(example_tr,4),"_(mo)",sep="",collapse="") # change from _alt to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[2]]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
            }else if(example_tr %in% tr$V1){ # if the name matches --> there is not alternative 
                
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- tr_line$V4[1]
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
                
                
            }else{
                example_tr <- paste(RemoveElements(example_tr,0),"_(mo)",sep="",collapse="") # change to _(mo)
                tr_line <- tr[tr$V1==example_tr,]
                exon_number_seq_i <- tr_line$V2[1]
                dna_seq_i <- tr_line$V3[1]
                orf_seq_i <- uncollapse(tr_line$V4[1])[[1]][[1]] # but now extrac the first ORF
                
                exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
                dna_seq <- c(dna_seq,dna_seq_i)
                orf_seq <- c(orf_seq,orf_seq_i)
            }
            
        }
        
        # add to comp_tr 
        
        comp_tr <- taRifx::remove.factors(as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq)))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        suppressWarnings(if(excluded_exons == 0){
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[,4])
                gene_end <- max(exon_data[,5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }else{
            if(mult==as.numeric(-1)){
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                
                
                
            }else{
                gene_start <- min(exon_data[-c(excl_exons_save),4])
                gene_end <- max(exon_data[-c(excl_exons_save),5])
                gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
                text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
                text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
                
            }
        }
        )
        
        
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            if(is.na(comp_tr$ORF[n])){
                #do nothing
            }else{
                ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # make x == n linenumber of this transcript
                
                x <- n
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                if(ORFS_only==TRUE){
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x],lty = 2)
                            
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[orf_exon_first]
                            max_right_ex <- max_left+exons_fitted_save$V5[orf_exon_last]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            
                            
                        }
                        
                        if(i %in% c(orf_exon_last,orf_exon_first)){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                            
                        }else if(i%in%UTR_exons){
                            #rect(xleft = max_left+exons$V4[i], xright = max_left+exons$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)
                            
                        }else{
                            
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        }
                        
                        
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                }else{
                    
                    for (i in ex_seq){  # loop over exon sequence
                        i <- as.numeric(i)
                        
                        if (i == min(as.numeric(ex_seq))){ # plot exon line
                            max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                            max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                            rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                            #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                            
                        }
                        
                        if(i%in%UTR_exons){
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                        }else{
                            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        }
                        
                        text(x = max_left+exons_fitted_save$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                        
                        if(orf_exon_first==orf_exon_last){
                            rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }else{
                            if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                                rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }
                            
                            if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                                rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                                
                            }}
                        
                        
                    }
                } #closing from if ORFS only else 
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }}
        
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        col_table_perc <- colorRampPalette(colors = c('grey',"darkred",'darkblue'))
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 10200+additional_place/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('darkblue','white','darkred'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = 9500, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
                
                for(i in allrows(comp_tr)){
                    # i <- 4
                    for(n in c(1:num_tissues)){
                        #n <- 1
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            if(val_transf==0){
                                perc_col <- col_table_perc(1000)[1]
                            }else{
                                perc_col <- col_table_perc(1000)[val_transf]
                            }
                            rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        }
                        
                        
                    }
                }
                # legend 
                rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-additional_place/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 10200+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
                
                
            }
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-7):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            if(FPKM == TRUE){
                text(labels = "Mean FPKM",x = 0+additional_place_2/2, y =ref[1]+200, cex = 1.5)
                
                col_table_fpkm <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        fpkm_val <- tissue_df[i,(2*n)]
                        if(is.na(fpkm_val)){
                            fpkm_col <- "white"
                        }else if(fpkm_val<=0.1){
                            fpkm_col <- col_table_fpkm(7)[1]
                        }else if(fpkm_val<=0.5){
                            fpkm_col <- col_table_fpkm(7)[2]
                        }else if(fpkm_val<=0.999){
                            fpkm_col <- col_table_fpkm(7)[3]
                        }else if(fpkm_val<=1.5){
                            fpkm_col <- col_table_fpkm(7)[4]
                        }else if(fpkm_val<=2.5){
                            fpkm_col <- col_table_fpkm(7)[5]
                        }else if(fpkm_val<=5){
                            fpkm_col <- col_table_fpkm(7)[6]
                        }else if(fpkm_val>5){
                            fpkm_col <- col_table_fpkm(7)[7]
                        }
                        
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = fpkm_col, border = TRUE)
                    }
                }
                # legend below
                legend(legend = c("Not found","<0.1","<0.5","<1","<1.5","<2.5","<5",">5"), x = -1300, y = 1450, fill=c("white",col_table_fpkm(7)[1:7]), title = "Heatmap legend", cex = 0.8) 
                
            }else{
                text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
                
                
                #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
                
                for(i in allrows(comp_tr)){
                    for(n in c(1:num_tissues)){
                        
                        perc_val <- tissue_df[i,(2*n)]
                        if(is.na(perc_val)){
                            perc_col <- "white"
                        }else{
                            val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                            perc_col <- col_table_perc(1000)[val_transf]
                            rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                            
                        }
                        
                    }
                }
                # legend below
                
                rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
                part<-1200/100
                for(i in c(1:100)){
                    # i <- 1
                    n <- i-1
                    rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
                }
                text(labels = "0%",x = 0+75, y =ref[2]+10, cex = 0.75)
                text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
                
                #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
                
                
            }
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# gapless
ExonCoveragePlot <-function(cov_file, distance = 100, interval = 100, log = FALSE){
    # cov file = output from samtools depth
    
    default_dist <- distance
    int_value <- interval
    
    cov_df <- read.table(file = cov_file, stringsAsFactors = FALSE)
    
    cov_df$V2 <- cov_df$V2 - min(cov_df$V2)
    
    for(i in allrows(cov_df)){
        n <- i+1
        if(n <= nrow(cov_df)){
            if((cov_df$V2[n]-cov_df$V2[i])>=default_dist){
                distance <- cov_df$V2[n]-cov_df$V2[i]
                subtr_dist <- distance-default_dist
                for(x in n:max(allrows(cov_df))){
                    cov_df$V2[x] <- cov_df$V2[x] - subtr_dist
                }
            }
        }
        
    }
    
    
    intervalls <- allrows(cov_df)[seq(1,length(allrows(cov_df)),int_value)]
    max_val <- length(allrows(cov_df))
    
    if(log == TRUE){
        cov_df$V3 <- log10(cov_df$V3+1)
    }
    
    
    if(int_value>=2){
        for(i in intervalls){
            #i <- intervalls[length(intervalls)]
            if(i == 1){
                n <- 1
                val <- mean(cov_df[i:(i+int_value),3])
                interval_df <- cbind(NA,NA)
                interval_df <- as.data.frame(rbind(interval_df,(c(n,val))))
                interval_df <- interval_df[-1,]
            }else if(i>=(max_val-(int_value-1))){
                n <- n+1
                val <- mean(cov_df[i:max_val,3])
                interval_df <- as.data.frame(rbind(interval_df,c(n,val)))
            }else{
                n <- n+1
                val <- mean(cov_df[i:(i+int_value),3])
                interval_df <- as.data.frame(rbind(interval_df,c(n,val)))
            }
        }
        
        plot(interval_df$V2 ~ interval_df$V1, type = "l", ylab = "coverage", xaxt='n', xlab = paste("resolution:",int_value,"bases", collapse=""), main = cov_file)
        
    }else{
        plot(cov_df$V3~cov_df$V2, type = "l",ylab = "coverage", xaxt='n', xlab = paste("resolution: 1 base"), main = cov_file)
    }
    
}
# added non-default option of using log coverage values

ReadGTFApplyCutoff <-function(file, unit = "FPKM", relative = FALSE, cutoff = 0){
    
    # load gtf file
    
    gtf <- read.table(file, stringsAsFactors = FALSE, sep = "\t")
    
    gtf_transcripts <- gtf[gtf$V3 == "transcript",]
    transcript_table <- c()
    
    for(i in allrows(gtf_transcripts)){
        row_tr <- uncollapse(gtf_transcripts[i,9])[[1]][c(4,6,8,10)]
        transcript_table <- rbind(transcript_table,row_tr)
    }
    
    colnames(transcript_table) <- c("name","cov","FPKM","TPM")
    for(i in allcols(transcript_table)){ # remove ";"
        for(n in allrows(transcript_table)){
            transcript_table[n,i] <- RemoveElements(transcript_table[n,i],1)
        }
    }
    
    transcript_table <- taRifx::remove.factors(as.data.frame(transcript_table))
    
    # get cutoff specifications
    
    tr_table <- transcript_table[,1]
    tr_table <- taRifx::remove.factors(as.data.frame(cbind(tr_table,transcript_table[,unit])))
    tr_table$V2 <- as.numeric(tr_table$V2)
    
    highest <- order(decreasing = TRUE, as.numeric(tr_table$V2))[1]
    
    if(relative==TRUE){
        tr_table$V2 <- tr_table$V2/tr_table$V2[highest]
    }
    
    
    # apply cutoff to variants
    
    tr_table <- tr_table[c(tr_table$V2>=cutoff),]
    allowed_variants <- tr_table$tr_table
    
    # go over gtf input and remove everything containing unallowed variants
    keep <- c()
    for(i in allrows(gtf)){
        for(n in allowed_variants){
            if(grepl(pattern = n,x = gtf$V9[i])==TRUE){
                keep <- c(keep,i)
            }    
        }
        
    }
    
    gtf_out <- gtf[keep,]
    return(gtf_out)
}

CutOffORFs <- function(compared_transcripts, cutoff = 0.01){
    
    #original idea: include and subset transcripts based on results
    
    #transcripts <- tr_l1
    #compared_transcripts <- tr_compared
    
    #tr <- transcripts
    comp_tr <- compared_transcripts
    
    ### remove every transcript that does not occur with [cutoff] fraction in at least 1 tissue
    
    # first remove ORFless and references
    
    comp_tr <- comp_tr[comp_tr$ORF!="NO ORF FOUND",]
    
    keep <- c()
    for(i in allrows(comp_tr)){
        if(grepl(pattern = "_ref",x = comp_tr$all[i])){
            #nothing
        }else{keep <- c(keep,i)}
    }
    
    comp_tr <- comp_tr[keep,]
    
    
    # now for making it smaller
    
    above <- c() # above coverage
    
    cols <- c()
    for(i in 1:length(colnames(comp_tr))){
        if(grepl(pattern="mean",x=colnames(comp_tr)[i])){
            cols <- c(cols,i)
        }
    }
    
    for(i in allrows(comp_tr)){
        if(comp_tr$mean_all[i]>=cutoff){
            above <- c(above,i)
        }else{
            for(n in cols){
                if(is.na(comp_tr[i,n])){
                    #do nothing
                }else{
                    if(comp_tr[i,n]>=cutoff){
                        above <- c(above,i)
                    }}
            }
            
            
        }
        
        
    }
    
    above <- unique(above)
    
    comp_tr_new <- comp_tr[above,]
    
    return(comp_tr_new)
}

CutOffORFsAndRemoveTranscripts <- function(compared_transcripts, cutoff = 0.01, transcript_table){
    
    #comp_tr_cutoff <- CutOffORFs(compared_transcripts = tr_compared, cutoff = 0.01)
    comp_tr_cutoff <- CutOffORFs(compared_transcripts, cutoff = cutoff)
    
    ### now remove every transcript of transcript table that doesnt belong to comp_tr_cutoff
    
    #transcript_table <- transcripts
    # extract all transcripts from comp_tr_cutoff
    tr_names <- c()
    for(i in allrows(comp_tr_cutoff)){
        tr_names <- c(tr_names,uncollapse(comp_tr_cutoff$all[i])[[1]])
    }
    
    tr_names <- unique(tr_names)
    
    for(i in 1:length(tr_names)){
        if(grepl(pattern = "alt", x = tr_names[i])){
            tr_names[i] <- RemoveElements(tr_names[i],amount = 4)
        }
    }
    
    tr_names <- unique(tr_names)
    
    keep <- c()
    
    for(i in allrows(transcript_table)){
        if(grepl(pattern = "(mo)", x = transcript_table$V1[i])){
            variant_name <- RemoveElements(transcript_table$V1[i], amount = 5)
            if(variant_name%in%tr_names){
                keep <- c(keep,i)
            }
        }else{
            variant_name <- transcript_table$V1[i]
            if(variant_name%in%tr_names){
                keep <- c(keep,i)
            }
        }
    }
    
    transcripts_new <- transcript_table[keep,]
    rownames(transcripts_new) <- 1:length(rownames(transcripts_new)) 
    return(transcripts_new)
}

BuildExonTable <- function(exons,transcripts){
    
    ex <- exons
    tr <- transcripts
    # comp_tr <- tr_compared
    #comp_tr  <- compared_transcripts
    
    # transform transcripts
    
    tr <- AddSamplesAndRelFPKMToTranscripts(tr)
    
    
    # make new dataframe
    
    colnames <- c("transcript", "exons", rownames(ex))
    
    tr_table <- rbind(colnames)
    colnames(tr_table) <- colnames
    tr_table <- tr_table[-1,]
    
    fraction_in_sample <- c()
    
    for(i in allrows(tr)){
        #i <- 1
        new_line <- c(tr$V1[i],tr$V2[i])
        
        ex_seq <- uncollapse(new_line[2])[[1]][-1]
        
        for(n in allrows(exons)){
            #n <- 4
            if(n %in% ex_seq){
                new_line <- c(new_line,"1")
            }else{
                new_line <- c(new_line," ")
            }
        }
        
        
        fraction <- as.numeric(RemoveElements(uncollapse(tr$V7[i])[[1]][8], amount = 1))
        
        fraction_in_sample <- c(fraction_in_sample,fraction)
        tr_table <- rbind(tr_table,new_line)
        
    }
    
    rownames(tr_table) <- c(1:nrow(tr_table))
    
    #tr_table <- as.data.frame(cbind(tr_table[,1:2],fraction_in_sample,tr_table[3:length(colnames(tr_table))]))
    
    tr_table_w_fraction <- taRifx::remove.factors(as.data.frame(cbind(tr_table[,c(1,2)],fraction_in_sample,tr_table[,c(3:length(colnames(tr_table)))])))
    
    # later: add seq, coding seq & orf
    #tr_table_w_fraction <- taRifx::remove.factors(as.data.frame(cbind(tr_table_w_fractions,)))
    
    return(tr_table_w_fraction)
    sum(as.numeric(tr_table_w_fraction$fraction_in_sample))
    
}

MergeDuplicateTranscripts <- function(transcripts, exons = exons){
    
    # idea of this function: merge all transcripts not just based on orf, but based on exons
    # however do ignore the very 5' end of the very 5' exons and vice versa
    
    # dont forget - stranded genes
    #transcripts
    
    for(i in allrows(transcripts)){
        #i <- 1
        print(transcripts$V1[i])
        trans_i <- transcripts$V1[i]
        ex_i <- transcripts$V2[i]         
        if(i == 1){
            uniq_table <- as.data.frame(rbind(c(transcripts$V2[i],transcripts$V1[i])))
            x <- sapply(uniq_table, is.factor)
            uniq_table[x] <- lapply(uniq_table[x], as.character)
            
        }else{ #i <- 2
            equal <- FALSE
            for(n in allrows(uniq_table)){
                # i  <- 2
                if(uniq_table[n,1]==trans_i){
                    equal <- TRUE
                    #uniq_table$V1[n] <- c(uniq_table$V1[n],ex_i)
                    uniq_table$V2[n] <- paste(collapse = " ", sep = " ", c(uniq_table$V2[n],trans_i))
                }
            }
            
            
            if(equal == FALSE){
                found <- FALSE
                n <- 0
                while(found == FALSE){
                    n <- n + 1
                    first_match <- FALSE
                    last_match <- FALSE
                    
                    first_i <- strsplit(as.character(ex_i), split = " ")[[1]]
                    if(first_i[1]==""){
                        first_i <- first_i[2]
                    }else{
                        first_i <- first_i[1]
                    }
                    
                    first_n <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                    if(first_n[1]==""){
                        first_n <- first_n[2]
                    }else{
                        first_n <- first_n[1]
                    }
                    
                    last_n <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                    last_n <- last_n[length(last_n)]
                    last_i <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                    last_i <- last_i[length(last_i)]
                    
                    if(exons$V5[as.numeric(first_i)]==exons$V5[as.numeric(first_n)]){
                        first_match <- TRUE
                    }
                    
                    if(exons$V4[as.numeric(last_i)]==exons$V4[as.numeric(last_n)]){
                        last_match <- TRUE
                    }
                    
                    if(first_match == TRUE && last_match == TRUE){
                        # check if rest is equal
                        print("match possible")
                        ex_i_string <- strsplit(as.character(ex_i), split = " ")[[1]]
                        ex_n <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                        
                        if(ex_i_string[1]==""){
                            ex_i_string <- ex_i_string[-c(1:2)]
                        }else{
                            ex_i_string <- ex_i_string[-c(1)]
                        }
                        
                        if(ex_n[1]==""){
                            ex_n <- ex_n[-c(1:2)]
                        }else{
                            ex_n <- ex_n[-c(1)]
                        }
                        
                        ex_n <- ex_n[-length(ex_n)]
                        ex_i_string <- ex_i_string[-length(ex_i_string)]
                        
                        if(paste(ex_i_string, sep = " ", collapse = " ") == paste(ex_n,sep = " ", collapse = " ")){#
                            uniq_table$V2[n] <- paste(collapse = " ", sep = " ", c(uniq_table$V2[n],trans_i))
                            
                            found <- TRUE
                        }else if(n == nrow(uniq_table)){
                            uniq_table <- as.data.frame(rbind(uniq_table,c(ex_i, trans_i)))
                            x <- sapply(uniq_table, is.factor)
                            uniq_table[x] <- lapply(uniq_table[x], as.character)
                            found <- TRUE 
                        }else{
                            # do nothing
                        }
                        
                        
                    }else if(n == nrow(uniq_table)){
                        uniq_table <- as.data.frame(rbind(uniq_table,c(ex_i, trans_i)))
                        x <- sapply(uniq_table, is.factor)
                        uniq_table[x] <- lapply(uniq_table[x], as.character)
                        found <- TRUE 
                    }
                    
                    
                }
                
                
                
            }
            
            
        }
    }
    
    return(uniq_table)
    
    
}

MergeDuplicateTranscripts <- function(transcripts, exons = exons){
    
    # idea of this function: merge all transcripts not just based on orf, but based on exons
    # however do ignore the very 5' end of the very 5' exons and vice versa
    
    # dont forget - stranded genes
    #transcripts
    
    for(i in allrows(transcripts)){
        #i <- 39
        print(transcripts$V1[i])
        trans_i <- transcripts$V1[i]
        ex_i <- transcripts$V2[i]         
        if(i == 1){
            uniq_table <- as.data.frame(rbind(c(transcripts$V2[i],transcripts$V1[i])))
            x <- sapply(uniq_table, is.factor)
            uniq_table[x] <- lapply(uniq_table[x], as.character)
            
        }else{ #i <- 2
            equal <- FALSE
            for(n in allrows(uniq_table)){
                # i  <- 2
                if(uniq_table[n,1]==trans_i){
                    equal <- TRUE
                    #uniq_table$V1[n] <- c(uniq_table$V1[n],ex_i)
                    uniq_table$V2[n] <- paste(collapse = " ", sep = " ", c(uniq_table$V2[n],trans_i))
                }
            }
            
            
            if(equal == FALSE){
                found <- FALSE
                n <- 0
                while(found == FALSE){
                    n <- n + 1
                    # n <- 29
                    first_match <- FALSE
                    last_match <- FALSE
                    
                    first_i <- strsplit(as.character(ex_i), split = " ")[[1]]
                    if(first_i[1]==""){
                        first_i <- first_i[2]
                    }else{
                        first_i <- first_i[1]
                    }
                    
                    first_n <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                    if(first_n[1]==""){
                        first_n <- first_n[2]
                    }else{
                        first_n <- first_n[1]
                    }
                    
                    last_n <- strsplit(as.character(ex_i), split = " ")[[1]]
                    last_n <- last_n[length(last_n)]
                    last_i <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                    last_i <- last_i[length(last_i)]
                    
                    if(exons$V5[as.numeric(first_i)]==exons$V5[as.numeric(first_n)]){
                        first_match <- TRUE
                    }
                    
                    if(exons$V4[as.numeric(last_i)]==exons$V4[as.numeric(last_n)]){
                        last_match <- TRUE
                    }
                    
                    if(first_match == TRUE && last_match == TRUE){
                        # check if rest is equal
                        print("match possible")
                        ex_i_string <- strsplit(as.character(ex_i), split = " ")[[1]]
                        ex_n <- strsplit(as.character(uniq_table[n,1]), split = " ")[[1]]
                        
                        if(ex_i_string[1]==""){
                            ex_i_string <- ex_i_string[-c(1:2)]
                        }else{
                            ex_i_string <- ex_i_string[-c(1)]
                        }
                        
                        if(ex_n[1]==""){
                            ex_n <- ex_n[-c(1:2)]
                        }else{
                            ex_n <- ex_n[-c(1)]
                        }
                        
                        ex_n <- ex_n[-length(ex_n)]
                        ex_i_string <- ex_i_string[-length(ex_i_string)]
                        
                        if(paste(ex_i_string, sep = " ", collapse = " ") == paste(ex_n,sep = " ", collapse = " ")){#
                            uniq_table$V2[n] <- paste(collapse = " ", sep = " ", c(uniq_table$V2[n],trans_i))
                            
                            found <- TRUE
                        }else if(n == nrow(uniq_table)){
                            uniq_table <- as.data.frame(rbind(uniq_table,c(ex_i, trans_i)))
                            x <- sapply(uniq_table, is.factor)
                            uniq_table[x] <- lapply(uniq_table[x], as.character)
                            found <- TRUE 
                        }else{
                            # do nothing
                        }
                        
                        
                    }else if(n == nrow(uniq_table)){
                        uniq_table <- as.data.frame(rbind(uniq_table,c(ex_i, trans_i)))
                        x <- sapply(uniq_table, is.factor)
                        uniq_table[x] <- lapply(uniq_table[x], as.character)
                        found <- TRUE 
                    }
                    
                    
                }
                
                
                
            }
            
            
        }
    }
    
    return(uniq_table)
    
    
}

RelativeTranscriptComparison <- function(unique_transcripts = unique_transcripts, transcripts = transcripts, genename_length = 6){
    
    #tr <- unique_transcripts
    tr <- transcripts
    tr_uniq <- unique_transcripts
    
    # get relative values
    tr_rel <- AddSamplesAndRelFPKMToTranscripts(tr_table = tr, genename_length = genename_length)
    
    # convers unique_transcripts
    freqs <- c()
    for(i in allrows(tr_uniq)){
        #i <- 1
        freqs <- c(freqs,length(uncollapse(tr_uniq$V2[i])[[1]]))
    }
    tr_uniq <- cbind(tr_uniq,freqs)
    tr_uniq <- tr_uniq[order(tr_uniq$freqs)[length(tr_uniq$freqs):1],-3]
    
    # now add coloumns to tr_uniq
    
    tr_uniq <- cbind(tr_uniq,NA,0)
    colnames(tr_uniq)[3] <- "all"
    colnames(tr_uniq)[4] <- "all_rel"
    
    
    # how much coloumns do you need?
    
    tissues <- as.character(unique(tr_rel$type))
    
    len_new_cols <- length(tissues) * 2
    
    for(n in tissues){
        tr_uniq <- cbind(tr_uniq,NA,0)
        colnames(tr_uniq)[c((ncol(tr_uniq)-1):ncol(tr_uniq))] <- c(n,paste(n,"_rel",sep="",collapse=""))
    }
    
    # fill it wiht values and tissue dependency
    
    for (i in allrows(tr_uniq)){
        # i <- 1
        samples <- uncollapse(tr_uniq$V2[i])[[1]]
        # add and in the end divide
        for(n in samples){
            #n <- samples[1]
            tr_tmp <- tr_rel[tr_rel$V1==n,]
            tissue_tmp <- tr_tmp$type[1]
            rel_tmp <- as.numeric(RemoveElements(char = uncollapse(tr_tmp$V7[1])[[1]][8],amount = 1))
            
            col_tmp <- which( colnames(tr_uniq)==tissue_tmp )
            tr_uniq[i,col_tmp] <- paste(tr_uniq[i,col_tmp],n)
            tr_uniq[i,(col_tmp+1)] <- tr_uniq[i,(col_tmp+1)]+rel_tmp
            
        }
        
        
    }
    
    indiv <- unique(tr_rel$affil)
    tr_uniq$all_rel <- as.numeric(tr_uniq$all_rel)
    
    for(i in allrows(tr_uniq)){
        # i <- 1
        tr_uniq$all_rel[i] <- sum(as.data.frame(tr_uniq[i,grepl(x = colnames(tr_uniq), pattern = "_rel")]))
        tr_uniq$all_rel[i] <- tr_uniq$all_rel[i]/length(indiv)
        
    }
    
    
    for(i in tissues){
        # i <- tissues[1]
        col_tmp <- which( colnames(tr_uniq)==i )+1
        
        indiv_tmp <- indiv[grepl(pattern = i, x = indiv)]
        replicates_tmp <- length(unique(indiv_tmp))
        
        tr_uniq[,col_tmp] <- tr_uniq[,col_tmp]/replicates_tmp
        
    }
    
    
    for(i in tissues){
        # i <- tissues[2]
        col_tmp <- which( colnames(tr_uniq)==i )
        for(n in allrows(tr_uniq)){
            # n <- 1
            if(is.na(tr_uniq[n,i])){
                # dont do anything
            }else{
                tr_uniq[n,i] <- RemoveElements(tr_uniq[n,i], amount = 3, start = "first")
            }
            
        }
    }
    
    tr_uniq_finished <- tr_uniq[,-3]
    
    return(tr_uniq_finished)
    
}

CutOffTranscripts<- function(compared_transcripts, cutoff = 0.01){
    
    #original idea: include and subset transcripts based on results
    
    #transcripts <- tr_l1
    #compared_transcripts <- tr_compared
    
    #tr <- transcripts
    comp_tr <- compared_transcripts
    
    ### remove every transcript that does not occur with [cutoff] fraction in at least 1 tissue
    
    # first remove ORFless and references
    
    #comp_tr <- comp_tr[comp_tr$ORF!="NO ORF FOUND",]
    
    keep <- c()
    for(i in allrows(comp_tr)){
        if(grepl(pattern = "_ref",x = comp_tr[i,2])){
            #nothing
        }else{keep <- c(keep,i)}
    }
    
    comp_tr <- comp_tr[keep,]
    
    
    # now for making it smaller
    
    above <- c() # above coverage
    
    cols <- c()
    for(i in 1:length(colnames(comp_tr))){
        if(grepl(pattern="rel",x=colnames(comp_tr)[i])){
            cols <- c(cols,i)
        }
    }
    
    for(i in allrows(comp_tr)){
        if(comp_tr[i,3]>=cutoff){
            above <- c(above,i)
        }else{
            for(n in cols){
                if(is.na(comp_tr[i,n])){
                    #do nothing
                }else{
                    if(comp_tr[i,n]>=cutoff){
                        above <- c(above,i)
                    }}
            }
            
            
        }
        
        
    }
    
    above <- unique(above)
    
    comp_tr_new <- comp_tr[above,]
    
    return(comp_tr_new)
}

CutOffTranscriptsAndRemove<- function(compared_transcripts, cutoff = 0.01, transcript_table){
    
    #comp_tr_cutoff <- CutOffORFs(compared_transcripts = tr_compared, cutoff = 0.01)
    comp_tr_cutoff <- CutOffTranscripts(compared_transcripts, cutoff = cutoff)
    # comp_tr_cutoff <- tr_relative_cutoff
    
    ### now remove every transcript of transcript table that doesnt belong to comp_tr_cutoff
    #
    #transcript_table <- transcripts
    # extract all transcripts from comp_tr_cutoff
    tr_names <- c()
    for(i in allrows(comp_tr_cutoff)){
        tr_names <- c(tr_names,uncollapse(comp_tr_cutoff[i,2])[[1]])
    }
    
    tr_names <- unique(tr_names)
    
    
    keep <- c()
    
    for(i in allrows(transcript_table)){
        variant_name <- transcript_table$V1[i]
        if(variant_name%in%tr_names){
            keep <- c(keep,i)
        }
    }
    
    transcripts_new <- transcript_table[keep,]
    rownames(transcripts_new) <- 1:length(rownames(transcripts_new)) 
    return(transcripts_new)
}

TranscriptsWithHeatmaps <- function(exon_data, transcript_list, compared_transcripts, vis.ref.file, domains = FALSE, IDs = FALSE, exon_text = TRUE){
    
    ### heatmap colors
    
    col_table_perc <- colorRampPalette(colors = c("deepskyblue4","darkred"))
    
    
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            
            ##############
            
            
            
            example_tr <- uncollapse(comp_tr[i,2])[[1]][1]
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }
        
        # add to comp_tr 
        
        comp_tr <- as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons_fitted_save$V4 * mult
        exons$V4 <- exons_fitted_save$V4 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            # n <- 1
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # are there multiple orfs?
            if(grepl(x = comp_tr$V2[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
                
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }else{ # if there is one
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][1]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                ################# # give the second orf!!!
                
                
                
                
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][2]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                
                
                
                
                
                
            }
        }
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                    
                    
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
            
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        
                    }
                    
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
            
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}

TranscriptsWithHeatmaps <- function(exon_data, transcript_list, compared_transcripts, vis.ref.file, domains = FALSE, IDs = FALSE, exon_text = TRUE){
    
    ### heatmap colors
    
    col_table_perc <- colorRampPalette(colors = c("deepskyblue4","darkred"))
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            
            ##############
            
            
            
            example_tr <- uncollapse(comp_tr[i,2])[[1]][1]
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }
        
        # add to comp_tr 
        
        comp_tr <- as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons$V4 * mult
        exons$V5 <- exons$V5 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            # n <- 4
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # are there multiple orfs?
            if(grepl(x = comp_tr$V2[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
                
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        
                        # here lies the problem
                        
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }else{ # if there is one
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][1]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                ################# # give the second orf!!!
                
                
                
                
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][2]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                
                
                
                
                
                
            }
        }
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                    
                    
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
            
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        
                    }
                    
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
            
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}

TranscriptsWithHeatmaps <- function(exon_data, transcript_list, compared_transcripts, vis.ref.file, domains = FALSE, IDs = FALSE, exon_text = TRUE){
    
    ### heatmap colors
    
    col_table_perc <- colorRampPalette(colors = c("deepskyblue4","darkred"))
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            
            ##############
            
            
            
            example_tr <- uncollapse(comp_tr[i,2])[[1]][1]
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }
        
        # add to comp_tr 
        
        comp_tr <- as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons$V4 * mult
        exons$V5 <- exons$V5 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            # n <- 4
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # are there multiple orfs?
            if(grepl(x = comp_tr$V2[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
                
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        
                        # here lies the problem
                        
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }else{ # if there is one
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][1]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                ################# # give the second orf!!!
                
                
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][2]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                
                
                
                
                
                
            }
        }
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                    
                    
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
            
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        
                    }
                    
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
            
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# now minus stranded genes work aswell: second orf is displayed correctly AND the first and last exon of the orf are displayed
TranscriptsWithHeatmaps <- function(exon_data, transcript_list, compared_transcripts, vis.ref.file, domains = FALSE, IDs = FALSE, exon_text = TRUE){
    
    ### heatmap colors
    
    col_table_perc <- colorRampPalette(colors = c("deepskyblue4","darkred"))
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            
            ##############
            
            
            
            example_tr <- uncollapse(comp_tr[i,2])[[1]][1]
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }
        
        # add to comp_tr 
        
        comp_tr <- as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons$V4 * mult
        exons$V5 <- exons$V5 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            # n <- 2
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # are there multiple orfs?
            if(comp_tr$orf_seq[n]=="NO ORF FOUND"){
                
                UTR_exons <- c(ex_seq) # define them all as utr exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    
                }
                
                
            }else if(grepl(x = comp_tr$V2[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
                
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        
                        # here lies the problem
                        
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }else{ # if there is one
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][1]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                ################# # give the second orf!!!
                
                
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][2]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                
                
                
                
                
                
            }
        }
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                    
                    
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
            
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        
                    }
                    
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
            
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# now accepts orf-less variants
TranscriptsWithHeatmaps <- function(exon_data, transcript_list, compared_transcripts, vis.ref.file, domains = FALSE, IDs = FALSE, exon_text = TRUE,palette=c("deepskyblue4","gold1","darkred")){
    
    ### heatmap colors
    
    col_table_perc <- colorRampPalette(colors = palette)
    
    ### load data & give new name
    
    exon_df <- exon_data
    tr <- transcript_list
    comp_tr <- compared_transcripts
    rownames(comp_tr) <- c(1:nrow(comp_tr))
    comp_tr_save <- comp_tr
    rownames(comp_tr) <- c(allrows(comp_tr))
    excluded_exons <- FALSE
    
    # does a domain_df exist?
    
    suppressWarnings(if(domains==FALSE){
        exist.domains<-"no"  
    }else{exist.domains<-"yes"})
    
    # extract gene_name
    tr_name <- uncollapse(tr$V1[1])[[1]]
    x <- "X"
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    tr_name <- RemoveElements(char = tr_name, amount = 1)
    x <- s2c(tr_name)[length(s2c(tr_name))]
    
    while(x!="_"){
        tr_name <- RemoveElements(char = tr_name, amount = 1)
        x <- s2c(tr_name)[length(s2c(tr_name))]
    }
    
    gene_name <- RemoveElements(char = tr_name, amount = 1)
    
    ### check if minus strand
    
    # check if it is a reverse-stranded gene and set mult = -1 (to multiply the parameters with -1 later)
    
    if(unique(exons$V7)=="-"){
        mult <- as.numeric(-1)
    }else if(unique(exons$V7)=="+"){
        mult <- as.numeric(1)}else{print("Different strands given!")}
    
    
    # figure out how much plots we need for showing every transcript in comp_tr
    
    num_rows <- allrows(comp_tr)
    set_size <- 9
    ratio <- max(num_rows)/set_size
    if((ratio-floor(ratio))==0){
        plot_replications <- floor(ratio)
    }else{
        plot_replications <- floor(ratio)+1
    }
    reps <- c()
    for(i in 1:plot_replications){
        # i <- 6
        if(i == plot_replications){
            tmp <- (set_size*(i-1)+1):(max(num_rows))
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }else{
            tmp <- (set_size*i-(set_size-1)):(set_size*i)
            reps <- c(reps,paste(tmp, collapse=" ",sep=" "))
        }
    }
    
    ### for each rep set of 9 transcripts
    
    for(i in reps){
        #i <- reps[1] 
        #print(i)
        rows_group <-uncollapse(i)[[1]]
        comp_tr <- comp_tr_save[rows_group,]    
        
        ### extract information about each compared transcript from transcript_list and add it to comp_tr
        
        exon_number_seq <- c()
        dna_seq <- c()
        orf_seq <- c()
        
        for(i in allrows(comp_tr)){
            # i <- 1      
            
            ##############
            
            
            
            example_tr <- uncollapse(comp_tr[i,2])[[1]][1]
            
            tr_line <- tr[tr$V1==example_tr,]
            exon_number_seq_i <- tr_line$V2[1]
            dna_seq_i <- tr_line$V3[1]
            orf_seq_i <- tr_line$V4[1]
            
            exon_number_seq <- c(exon_number_seq,exon_number_seq_i)
            dna_seq <- c(dna_seq,dna_seq_i)
            orf_seq <- c(orf_seq,orf_seq_i)
            
            
        }
        
        # add to comp_tr 
        
        comp_tr <- as.data.frame(cbind(comp_tr,exon_number_seq,dna_seq,orf_seq))
        
        ### exon preprocessing
        
        # load vis ref file
        exons <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        exons_fitted_save <- read.table(file = vis.ref.file,  sep = "\t",stringsAsFactors = FALSE )
        
        exons_fitted_save$V4 <- exons_fitted_save$V4 * mult
        exons_fitted_save$V5 <- exons_fitted_save$V5 * mult
        
        exons$V4 <- exons$V4 * mult
        exons$V5 <- exons$V5 * mult
        
        
        ### plot: basic parameters (don't forget minus strand)
        
        ylim_num <- 1300 # highest point of plot
        ref <- c(50+900,100+900) # first line
        addtoref <- c(-100,-200,-300,-400,-500,-600,-700,-800,-900) # next 9 lines
        max_left <- 100 * mult # leftmost coordinate
        max_right <- 10100 * mult# rightmost coordinate
        col_ref <- "darkblue"    
        
        if(mult==as.numeric(-1)){
            additional_place <- 200 # important for heatmap - develop that for - strand genes later!
            additional_place_2 <- 1200
        }else{additional_place <- 1200
        additional_place_2 <- 0}
        
        plot_x1 <- (additional_place+10200) * mult
        plot_x2 <- (-300-additional_place_2) * mult
        plot(x = c(plot_x1,plot_x2), y = c(0,0), ylim = c(-100,ylim_num+100), xlab ="",axes = FALSE, ylab = "", type = "l", col = "gray")
        
        
        
        ### plot: DNA locus and overlap exons (don't forget to use excluded dataframe)
        
        rect(xleft = max_left, xright = max_right, ybot = ref[1]+25, ytop = ref[2]-25) # plot first intron line
        
        # plot every possible exon in the highest line
        
        for (i in allrows(exons_fitted_save)){
            rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1], ytop = ref[2], col = col_ref,  border = FALSE)
        }
        
        #   add exonic composition of the genetic region 
        
        ex_sorted <- as.data.frame(exons_fitted_save)
        
        addtorefexon <-c(10,27,44,61,78,95,112,129)
        
        keep <- allrows(ex_sorted)
        x <- 1
        while(length(keep)>=1){
            new_keep <- c()
            for(i in keep){
                if(i == keep[1]){
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    end <- abs(ex_sorted$V5[i])
                }else if((abs(ex_sorted$V4[i])-end)>=15){
                    #text(x = max_left+exons$V4[i]+20, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.65)
                    rect(xleft = max_left + ex_sorted$V4[i], xright = max_left + ex_sorted$V5[i], ybot = ref[2]+addtorefexon[x], ytop = ref[2]+addtorefexon[x]+10, col = "darkblue", border = FALSE)
                    end <- abs(ex_sorted$V5[i]) 
                }else{new_keep <- c(new_keep,i)}
                
            }
            end <- c()
            keep <- new_keep
            x <- x +1
        }
        
        # add descriptions
        
        if(mult==as.numeric(-1)){
            text(x = max_right - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_right - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            text(x = max_right - 450, y =  ref[1]-15, labels = "(reverse strand)", col = col_ref, font = 3, cex = 0.75)
            
        }else{
            text(x = max_left - 450, y =  ref[1]+25, labels = "DNA locus", col = col_ref, font = 2)
            text(x = max_left - 450, y =  ref[2]+25, labels = "exons", col = "DARKBLUE", font = 3)
            
        }
        
        # plot: domains
        suppressWarnings(if(domains == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
        }else{
            
            dom <- domains
            
            for(i in allrows(dom)){
                ex_seq <- as.numeric(uncollapse(as.character(dom[i,2]))[[1]])
                
                
                firstexon <- floor(ex_seq[1]) # get integer exon
                perc_f_exon <- ex_seq[1]-firstexon# get leftmost coordinate
                left_coord <- exons_fitted_save$V4[firstexon]+((exons_fitted_save$V5[firstexon]-exons_fitted_save$V4[firstexon])*perc_f_exon)
                
                
                lastexon <- floor(ex_seq[length(ex_seq)]) # get integer exon
                perc_l_exon <- ex_seq[length(ex_seq)]-lastexon# get rightmost coordinate
                if(ex_seq[length(ex_seq)]%%1 > 0){
                    right_coord <- exons_fitted_save$V4[lastexon]+((exons_fitted_save$V5[lastexon]-exons_fitted_save$V4[lastexon])*perc_l_exon)
                }else{right_coord <- exons_fitted_save$V5[lastexon]}
                
                ex_seq[1] <- firstexon
                ex_seq[length(ex_seq)] <- lastexon
                if(length(unique(ex_seq))==1){
                    rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }else{for(n in ex_seq){
                    if(n==firstexon){
                        rect(xleft = max_left+left_coord, xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else if(n==lastexon){
                        rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }else{rect(xleft = max_left+exons_fitted_save$V4[n], xright = max_left+exons_fitted_save$V5[n], ybot = ref[1]-1000, ytop = ref[1], col = "grey", border = FALSE)# remove the black line after
                    }
                    
                    mean_dom <- max_left + (left_coord+right_coord)/2
                    text(labels = (paste(dom[i,1],sep="")), x = mean_dom, y =ref[2]-1095, cex = 0.8, srt = 45 )
                    
                }}
                
                
                rect(xleft = max_left+left_coord, xright = max_left+right_coord, ybot = ref[1]-1000, ytop = ref[1]-1000)
                
            }
        })
        
        
        
        
        # add absolute location
        
        chr <- exon_data$V1[1]
        
        if(mult==as.numeric(-1)){
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_right + 300, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            text(x = max_left - 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            
            
            
        }else{
            gene_start <- min(exon_data[,4])
            gene_end <- max(exon_data[,5])
            gene_start <- paste(chr,format(gene_start,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            gene_end <- paste(chr,format(gene_end,big.mark=",",scientific=FALSE),sep=":",collapse=":")
            text(x = max_left + 250, y =  ref[1]-20, labels = gene_start, col = col_ref, font = 2, cex = 0.90)
            text(x = max_right - 250, y =  ref[1]-20, labels = gene_end, col = col_ref, font = 2, cex = 0.90)
            
        }
        
        
        # text(labels = "Selection of splice variants | Source: RNA-Seq data of murine adipocytes", x = 1750, y = -30)
        
        text(labels = (paste(gene_name,' - splice variants',sep="")), x = 5100*mult, y = ylim_num -50, cex = 2)
        
        ### plot: comp_tr
        
        # basic color
        
        col_cov_utr <- "darkkhaki"
        col_cov_orf <- "darkolivegreen4"
        
        # careful with first & last exon
        # n <- 4
        
        for (n in allrows(comp_tr)){  # for chosen transcript rows
            # n <- 2
            ex_seq <- uncollapse(comp_tr$exon_number_seq[n])[[1]][-1]  #get exon seuqnece
            if(mult == as.numeric(-1)){ # if it is a reverse strand gene
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            # make x == n linenumber of this transcript
            
            x <- n
            
            # are there multiple orfs?
            if(comp_tr$orf_seq[n]=="NO ORF FOUND"){
                
                UTR_exons <- c(ex_seq) # define them all as utr exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    
                }
                
                
            }else if(grepl(x = comp_tr$V2[n], pattern = "(mo)")==FALSE){ # if there is NO alternative ORF
                
                
                # get exons which are in the ORF and in the UTR
                
                sequence <- as.character(comp_tr$dna_seq[n])
                orf <- as.character(comp_tr$orf_seq[n])
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                
                
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        
                        # here lies the problem
                        
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                #if(mult==as.numeric(-1)){
                #    text(x = max_right - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n])    
                #}else{
                #    text(x = max_left - 450, y =  ref[1]+addtoref[x]+25, labels = tr$V1[n]) 
                #}
            }else{ # if there is one
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][1]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                ################# # give the second orf!!!
                
                
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                orf <- as.character(uncollapse(comp_tr$orf_seq[n])[[1]][2]) # give him the first (longest)
                sequence <- as.character(comp_tr$dna_seq[n])    
                
                orf_position <- matchPattern(orf, sequence) # find orf in sequence
                
                orf_start_rel <- start(orf_position)/length(s2c(sequence)) # get relative orf start
                orf_end_rel <- end(orf_position)/length(s2c(sequence)) # get relative orf end
                
                exon_sum <- 0
                for (i in ex_seq){ # get exon sum 
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                orf_start <- orf_start_rel * exon_sum
                orf_end <- orf_end_rel * exon_sum
                
                exon_sum_max <- exon_sum
                exon_sum <- 0
                orf_exon_first <- NA
                orf_exon_last <- NA
                UTR_exons <- c()
                
                for(i in ex_seq){ # find exon in which the orf begins
                    #i <- 190
                    i <- as.numeric(i)
                    if(is.na(orf_exon_first)){
                        exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                        if(exon_sum<=orf_start){
                            UTR_exons <- c(UTR_exons,i)
                        }else if(is.na(orf_exon_first)){
                            orf_exon_first <- i
                        }}
                }
                
                ex_seq_back <- ex_seq[length(ex_seq):1] # make ex_seq backwards for finding orf end
                
                exon_sum <- exon_sum_max
                for(i in ex_seq_back){
                    i <- as.numeric(i)
                    exon_sum <- exon_sum - (abs(exons$V5[i])-abs(exons$V4[i]))
                    if(exon_sum>=orf_end){
                        UTR_exons <- c(UTR_exons,i)
                    }else if(is.na(orf_exon_last)){
                        orf_exon_last <- i
                    }
                }
                
                
                # group all exons which are not (or not completely) part of the orf  
                
                UTR_exons <- c(orf_exon_first, orf_exon_last, UTR_exons)
                
                # special case: the two exons in which the orf starts/end - find position within exon
                
                # orf start
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_first)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before first (incl. first)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                }
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_first])-abs(exons$V4[orf_exon_first])) # get sum of exons before first
                orf_start_exon <- orf_start-sum_before # subtract this from the orf start
                # now we have the value within the first exon
                
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_start_exon <-  orf_start_exon + exons$V5[orf_exon_first]
                    orf_start_exon_end <- exons$V4[orf_exon_first]      
                }else{
                    orf_start_exon <- orf_start_exon + exons$V4[orf_exon_first] # start for orf-start exon
                    orf_start_exon_end <- exons$V5[orf_exon_first]      # end for orf-start exon
                }
                
                # now for orf end
                
                ex_seq_tmp <- ex_seq[1:which(ex_seq == orf_exon_last)]
                
                exon_sum <- 0
                for(i in ex_seq_tmp){ # get exon sum of all exons before after last (incl. last)
                    i <- as.numeric(i)
                    exon_sum <- exon_sum + (abs(exons$V5[i])-abs(exons$V4[i]))
                } 
                
                sum_before <- exon_sum - (abs(exons$V5[orf_exon_last])-abs(exons$V4[orf_exon_last])) # get sum of exons after last
                orf_end_exon <- orf_end-sum_before # subtract this from the orf end
                # now we have the value within the last exon
                
                if(mult==as.numeric(-1)){ # if it is a -strand gene, you have to reverse v4 & v5
                    orf_end_exon <- orf_end_exon + exons$V5[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V5[orf_exon_last]      # start for orf-end exon      
                }else{
                    orf_end_exon <- orf_end_exon + exons$V4[orf_exon_last] # end for orf-end exon
                    orf_end_exon_start<- exons$V4[orf_exon_last]      # start for orf-end exon
                }
                
                
                # plot exons
                if(mult == as.numeric(-1)){ # if it is a reverse strand gene --> reverse back
                    ex_seq <- ex_seq[length(ex_seq):1]
                }
                
                # i <- 73
                for (i in ex_seq){  # loop over exon sequence
                    i <- as.numeric(i)
                    
                    if (i == min(as.numeric(ex_seq))){ # plot exon line
                        max_left_ex <- 10+max_left+exons_fitted_save$V4[min(as.numeric(ex_seq))]
                        max_right_ex <- max_left+exons_fitted_save$V5[max(as.numeric(ex_seq))]-10
                        #rect(xleft = max_left_ex, xright = max_right_ex, ybot = ref[1]+25+addtoref[x], ytop = ref[2]-25+addtoref[x])
                        #rect(xleft = max_left-10, xright =max_left+exons$V4[i], ybot = ref[1]+20+addtoref[x], ytop = ref[2]-20+addtoref[x], col = "white", border = NA)# remove the black line after
                        
                    }
                    
                    if(i%in%UTR_exons){
                        #rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+15+addtoref[x], ytop = ref[2]-15+addtoref[x], col = col_cov_utr, border = FALSE)    
                    }else{
                        rect(xleft = max_left+exons_fitted_save$V4[i], xright = max_left+exons_fitted_save$V5[i],ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                    }
                    
                    if(exon_text==TRUE){
                        text(x = max_left+exons_fitted_save$V4[i]+5, y = ref[1]+addtoref[x]-20, labels = i, cex = 0.45)
                    }
                    
                    if(orf_exon_first==orf_exon_last){
                        rect(xleft = max_left+orf_start_exon, max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                        
                    }else{
                        if(i == orf_exon_first){ # if it is the orf start, make additional "thick" orf start
                            rect(xleft = max_left+orf_start_exon, xright = max_left+orf_start_exon_end,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }
                        
                        if(i == orf_exon_last){ # if it is the orf end, make additional "thick" orf end
                            rect(xleft = max_left+orf_end_exon_start, xright = max_left+orf_end_exon,ybot = ref[1]+addtoref[x], ytop = ref[2]+addtoref[x], col = col_cov_orf, border = FALSE)
                            
                        }}
                    
                    
                }
                
                
                
                
                
                
                
            }
        }
        ### add IDs for each transcript
        
        suppressWarnings(if(IDs == FALSE){
            # text(labels = "Selection of splice variants | Source: RNA-Seq data", x = 1750, y = -30)
            
            #IDs <- c("a","b","","d","e","f","g","")
        }else{
            if(mult==as.numeric(-1)){
                for(i in c(1:length(IDs))){
                    
                    text(x = max_right - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
            }else{
                for(i in c(1:length(IDs))){
                    text(x = max_left - 450, y =  ref[1]+addtoref[i]+25, labels = IDs[i], col = "black", font = 1)
                }
                
            }
            
        })
        
        
        
        
        ### plot: heatmap 
        
        if(mult==as.numeric(1)){
            # define the room
            
            #rect(xleft = 10200, xright = 10200+additional_place, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 10200+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            
            text(labels = "Mean fraction",x = 10200+additional_place/2, y =ref[1]+125, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('darkblue',"white",'darkred'))
            
            for(i in allrows(comp_tr)){
                # i <- 4
                for(n in c(1:num_tissues)){
                    #n <- 1
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        if(val_transf==0){
                            perc_col <- col_table_perc(1000)[1]
                        }else{
                            perc_col <- col_table_perc(1000)[val_transf]
                        }
                        rect(xright = 10200+(n*tissue_place), xleft = 10200+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                    }
                    
                    
                }
            }
            # legend 
            rect(xright = 10200, xleft = 10200+additional_place, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-additional_place/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=10200+(n*part), xleft = 10200+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 10200+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 10200+additional_place-95, y =ref[2]+10, cex = 0.75)
            
            
            
        }else if(mult==as.numeric(-1)){ # if mult == -1
            
            # define the room
            
            # rect(xleft = 0, xright = 0+additional_place_2, ybottom = ref+addtoref[8], ytop = ref)
            
            # divide it by the number of differents tissues
            
            constant_cols <- c(1:3,(length(colnames(comp_tr))-2):length(colnames(comp_tr)))
            tissue_df <- comp_tr[,-c(constant_cols)] # now only the tissues and quantifications are remaining
            
            
            
            tissue_cols <-seq(1,max(allcols(tissue_df)),by=2)
            num_tissues <- length(tissue_cols)
            tissues <- colnames(tissue_df)[tissue_cols]
            tissue_place <- additional_place_2/num_tissues
            
            for(x in tissue_cols){
                for(y in allrows(tissue_df)){
                    if(is.na(tissue_df[y,x])){
                        tissue_df[y,x+1] <- NA
                    }
                }
            }
            
            
            
            for(i in c(1:num_tissues)){
                text(labels = tissues[i], x = 0+(i*tissue_place)-200, y =ref[2]-1050, cex = 1.2, srt = 45 )
            }
            
            text(labels = "Mean fraction",x = 0+additional_place_2/2, y =ref[1]+150, cex = 1.5)
            
            
            #col_table_perc <- colorRampPalette(colors = c('black',"darkorange",'darkgreen'))
            
            for(i in allrows(comp_tr)){
                for(n in c(1:num_tissues)){
                    
                    perc_val <- tissue_df[i,(2*n)]
                    if(is.na(perc_val)){
                        perc_col <- "white"
                    }else{
                        val_transf <- round(as.numeric(perc_val)*1000,digits = 0)
                        perc_col <- col_table_perc(1000)[val_transf]
                        rect(xright = 0+(n*tissue_place), xleft = 0+((n-1)*tissue_place), ybot = ref[2]+addtoref[i], ytop = ref[1]+addtoref[i], col = perc_col, border = TRUE)
                        
                    }
                    
                }
            }
            # legend below
            
            rect(xright = 0, xleft = 1200, ybot = ref[1]+10, ytop = ref[2]-10, col = "white", border = TRUE)
            part<-1200/100
            for(i in c(1:100)){
                # i <- 1
                n <- i-1
                rect(xright=0+(n*part), xleft = 0+(i*part) , ybot = ref[1]+10, ytop = ref[2]-10 , col = col_table_perc(100)[i], border = FALSE)
            }
            text(labels = "1%",x = 0+75, y =ref[2]+10, cex = 0.75)
            text(labels = "100%",x = 1200-100, y =ref[2]+10, cex = 0.75)
            
            #legend(legend = c("Not found","0-0.1","0.4-0.5","0.9-1"), x = -1300, y = 1450, fill=c("white",col_table_perc(10)[c(1,5,10)]), title = "Heatmap legend", cex = 0.8) 
            
            
            
            
            
        }
        
        ### legend
        if(mult==as.numeric(-1)){
            legend(legend = c("ORF","UTR"), x = -9000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
        }else{
            legend(legend = c("ORF","UTR"), x = 2000, y = ylim_num+120, fill=c(col_cov_orf, col_cov_utr), title = "Exon colours") 
            
        }
        
    } # loop over reps closing bracket
    
}
# new feature: palette colors as input
GetExonSequence <- function(references, exons, transcripts, locus){
    
    ### idea: find exons (and their locuses) of a transcript variants
    
    ref <- references
    x <- sapply(ref, is.factor)
    ref[x] <- lapply(ref[x], as.character)  
    # format references new
    ref_lines <- c(1:nrow(ref))
    ref_names <- c(ref_lines[ref_lines%%2==1])
    ref_tr <- c(ref_lines[ref_lines%%2==0])
    ref <- rbind(ref$V1[ref_names],ref$V1[ref_tr])
    
    ref <- as.data.frame(ref)
    x <- sapply(ref, is.factor)
    ref[x] <- lapply(ref[x], as.character)        
    
    # remove alignment gaps
    for(i in allcols(ref)){
        # i <- 2
        ref[2,i] <- gsub("-","",ref[2,i])
        
    }
    
    #gsub("-","","A-AADKT")
    
    
    
    ex <- exons
    
    tr <- transcripts
    
    loc <- locus
    
    # get the sequence of the exons
    
    ex_seqs <- c()
    
    for(i in allrows(exons)){
        #print(i)
        # i <- 1
        coords <- c(exons[i,c(4)],exons[i,c(5)])
        
        coords <- coords-exons[1,4]+1
        
        
        loc_i <- s2c(loc$V1[2])
        ex_seqs <- c(ex_seqs,c2s(loc_i[coords[1]:coords[2]]))
        
        
    }
    
    ex <- cbind(ex,ex_seqs)
    x <- sapply(ex, is.factor)
    ex[x] <- lapply(ex[x], as.character)        
    
    
    # search in transcripts wether ref is there --> then check the end
    equals <- c()
    for(i in allcols(ref)){
        ref_seq <-ref[2,i]
        # check if equal
        x <- FALSE
        for(n in allrows(transcripts)){
            if(transcripts$V3[n]==ref_seq){
                x <- transcripts$V1[n]     
            }
        }
        
        if(x == FALSE){
            equals <- c(equals,NA)
        }else{
            equals <- c(equals,x)
        }
        
    }
    
    
    within <- c()
    
    for(i in allcols(ref)){
        ref_seq <-ref[2,i]
        # check if equal
        x <- FALSE
        for(n in allrows(transcripts)){
            
            match_res <- matchPattern(pattern = tolower(ref_seq),subject = transcripts$V3[n])
            if(length(start(match_res))>=1){ # if there are matching patterns
                
                x <- transcripts$V1[n]          
            }
            
            
            
        }
        
        if(x == FALSE){
            within <- c(within,NA)
        }else{
            within <- c(within,x)
        }
        
    }
    
    
    # add equals and withins 
    
    
    ref <- rbind(ref[c(1,2),],equals,within)
    
    # for all withins --> get exons belonging
    
    for(i in allcols(ref)){
        print(i)
        if(!is.na(ref[4,i])){ # if there is one within transcript
            # i <- 3
            transcript_info <- transcripts[transcripts$V1==ref[4,i],c(2,3)]
            ex_seq <- transcript_info[1,1]
            tr_seq <- transcript_info[1,2]
            
            match_res <- matchPattern(pattern = tolower(ref[2,i]),subject = tr_seq)
            start <- start(match_res)
            end <- end(match_res)
            total <- length(s2c(tr_seq))
            
            ex_seq <- uncollapse(ex_seq)[[1]]
            if(ex_seq[1]==""){
                ex_seq <- ex_seq[-1]
            }
            
            # get exons that dont belong
            n <- 0
            first_found <- FALSE
            sum <- 0
            while(first_found == FALSE){
                n <- n+1
                ex_n <- as.numeric(ex_seq[n])
                sum <- sum  + (exons$V5[ex_n]-exons$V4[ex_n])+1
                if(sum >= start){
                    first_found <- ex_n
                    
                }
            }
            
            
            n <- 0
            last_found <- FALSE
            sum <- 0
            while(last_found == FALSE){
                length(s2c(tr_seq))
                n <- n+1
                ex_n <- as.numeric(ex_seq[n])
                sum <- sum  + (exons$V5[ex_n]-exons$V4[ex_n])+1
                if(sum >= end){
                    last_found <- ex_n
                    
                }
            }
            
            
            # add it to ref_table
            pos <- match(c(first_found,last_found),table = ex_seq)
            exon_string <- ex_seq[pos[1]:pos[2]]
            exon_string[1] <- paste("(",exon_string[1],")",sep = "", collapse = "")
            exon_string[length(exon_string)] <- paste("(",exon_string[length(exon_string)],")",sep = "", collapse = "")
            
            exon_string <- paste(exon_string, sep = " ", collapse = " ")
            ref[4,i] <- paste(ref[4,i],exon_string,sep = " ||| ", collapse = "||")
            
            
            
        }
    }
    
    
    
    
    # then look whicch exons exist within
    included_exons <- c()
    for(i in allcols(ref)){
        
        if(is.na(ref[4,i])){
            
            tr_seq <- ref[2,i]
            including <- c()
            for(n in allrows(ex)){
                #print(n)
                n_seq <- ex$ex_seqs[n]
                match_res <- matchPattern(pattern = n_seq,subject = tr_seq)
                start <- start(match_res)
                if(length(start)>0){
                    including <- c(including,n)
                }
                
                
            }
            
            incl_df <- c()
            
            for(n in including){
                incl_df <- rbind(incl_df,c(n,ex$V4[n],ex$V5[n]))
            }
            
            incl_df <- as.data.frame(incl_df)
            
            included_in_others <- c()
            for(n in allrows(incl_df)){
                # n <- 1
                num <- incl_df[n,1]
                reg_n <- c(incl_df[n,2],incl_df[n,3])
                group_table <- incl_df[incl_df$V2>=reg_n[1]&incl_df$V3<=reg_n[2]&incl_df$V1!=num,]  
                included_in_others <- c(included_in_others,group_table$V1)
                
                
            }
            
            included_in_others <- unique(included_in_others)
            
            unique_including <- including[-match(included_in_others,including)] 
            
            included_exons <- c(included_exons,paste(unique_including,sep=" ", collapse = " "))
        }else{
            included_exons <- c(included_exons,NA)
        }
        
    }
    
    
    ref <- as.data.frame(rbind(ref[c(1:4),],included_exons))
    
    
    return(ref)
    
}

GetExonSequence <- function(references, exons, transcripts, locus){
    
    print("may not work with minus stranded genes --> see comparison of refernces in adgrl2 in FIL")
    ### strand info
    
    if(unique(exons$V7)=="-"){
        strand <- "neg"
    }else if(unique(exons$V7=="+")){
        strand <- "pos"
    }
    
    
    ### idea: find exons (and their locuses) of a transcript variants
    
    ref <- references
    x <- sapply(ref, is.factor)
    ref[x] <- lapply(ref[x], as.character)  
    # format references new
    ref_lines <- c(1:nrow(ref))
    ref_names <- c(ref_lines[ref_lines%%2==1])
    ref_tr <- c(ref_lines[ref_lines%%2==0])
    ref <- rbind(ref$V1[ref_names],ref$V1[ref_tr])
    
    ref <- as.data.frame(ref)
    x <- sapply(ref, is.factor)
    ref[x] <- lapply(ref[x], as.character)        
    
    # remove alignment gaps
    for(i in allcols(ref)){
        # i <- 2
        ref[2,i] <- gsub("-","",ref[2,i])
        
    }
    
    ex <- exons
    
    tr <- transcripts
    
    loc <- locus
    
    # get the sequence of the exons
    
    ex_seqs <- c()
    
    for(i in allrows(exons)){
        #print(i)
        # i <- 1
        coords <- c(exons[i,c(4)],exons[i,c(5)])
        
        coords <- coords-exons[1,4]+1
        
        
        loc_i <- s2c(loc$V1[2])
        ex_seqs <- c(ex_seqs,c2s(loc_i[coords[1]:coords[2]]))
        
        
    }
    
    ex <- cbind(ex,ex_seqs)
    x <- sapply(ex, is.factor)
    ex[x] <- lapply(ex[x], as.character)        
    
    
    # search in transcripts wether ref is there --> then check the end
    equals <- c()
    for(i in allcols(ref)){
        ref_seq <-ref[2,i]
        # check if equal
        x <- FALSE
        for(n in allrows(transcripts)){
            if(transcripts$V3[n]==ref_seq){
                x <- transcripts$V1[n]     
            }
        }
        
        if(x == FALSE){
            equals <- c(equals,NA)
        }else{
            equals <- c(equals,x)
        }
        
    }
    
    
    within <- c()
    
    for(i in allcols(ref)){
        ref_seq <-ref[2,i]
        # check if equal
        x <- FALSE
        for(n in allrows(transcripts)){
            
            match_res <- matchPattern(pattern = tolower(ref_seq),subject = transcripts$V3[n])
            if(length(start(match_res))>=1){ # if there are matching patterns
                
                x <- transcripts$V1[n]          
            }
            
            
            
        }
        
        if(x == FALSE){
            within <- c(within,NA)
        }else{
            within <- c(within,x)
        }
        
    }
    
    
    # add equals and withins 
    
    
    ref <- rbind(ref[c(1,2),],equals,within)
    
    # for all withins --> get exons belonging
    
    for(i in allcols(ref)){
        print(i)
        # i <- 2
        if(!is.na(ref[4,i])){ # if there is one within transcript
            # i <- 3
            transcript_info <- transcripts[transcripts$V1==ref[4,i],c(2,3)]
            ex_seq <- transcript_info[1,1]
            tr_seq <- transcript_info[1,2]
            
            
            match_res <- matchPattern(pattern = tolower(ref[2,i]),subject = tr_seq)
            start <- start(match_res)
            end <- end(match_res)
            total <- length(s2c(tr_seq))
            
            ex_seq <- uncollapse(ex_seq)[[1]]
            
            if(ex_seq[1]==""){
                ex_seq <- ex_seq[-1]
            }
            
            if(strand=="neg"){
                ex_seq <- ex_seq[length(ex_seq):1]
            }
            
            
            # get exons that dont belong
            n <- 0
            first_found <- FALSE
            sum <- 0
            while(first_found == FALSE){
                n <- n+1
                ex_n <- as.numeric(ex_seq[n])
                sum <- sum  + (exons$V5[ex_n]-exons$V4[ex_n])+1
                if(sum >= start){
                    first_found <- ex_n
                    
                }
            }
            
            
            n <- 0
            last_found <- FALSE
            sum <- 0
            while(last_found == FALSE){
                length(s2c(tr_seq))
                n <- n+1
                ex_n <- as.numeric(ex_seq[n])
                sum <- sum  + (exons$V5[ex_n]-exons$V4[ex_n])+1
                if(sum >= end){
                    last_found <- ex_n
                    
                }
            }
            
            
            # add it to ref_table
            pos <- match(c(first_found,last_found),table = ex_seq)
            exon_string <- ex_seq[pos[1]:pos[2]]
            exon_string[1] <- paste("(",exon_string[1],")",sep = "", collapse = "")
            exon_string[length(exon_string)] <- paste("(",exon_string[length(exon_string)],")",sep = "", collapse = "")
            
            exon_string <- paste(exon_string, sep = " ", collapse = " ")
            ref[4,i] <- paste(ref[4,i],exon_string,sep = " ||| ", collapse = "||")
            
            
            
        }
    }
    
    
    
    
    # then look whicch exons exist within
    included_exons <- c()
    for(i in allcols(ref)){
        
        if(is.na(ref[4,i])){
            
            tr_seq <- ref[2,i]
            including <- c()
            for(n in allrows(ex)){
                #print(n)
                n_seq <- ex$ex_seqs[n]
                match_res <- matchPattern(pattern = n_seq,subject = tr_seq)
                start <- start(match_res)
                if(length(start)>0){
                    including <- c(including,n)
                }
                
                
            }
            
            incl_df <- c()
            
            for(n in including){
                incl_df <- rbind(incl_df,c(n,ex$V4[n],ex$V5[n]))
            }
            
            incl_df <- as.data.frame(incl_df)
            
            included_in_others <- c()
            for(n in allrows(incl_df)){
                # n <- 1
                num <- incl_df[n,1]
                reg_n <- c(incl_df[n,2],incl_df[n,3])
                group_table <- incl_df[incl_df$V2>=reg_n[1]&incl_df$V3<=reg_n[2]&incl_df$V1!=num,]  
                included_in_others <- c(included_in_others,group_table$V1)
                
                
            }
            
            included_in_others <- unique(included_in_others)
            
            unique_including <- including[-match(included_in_others,including)] 
            
            included_exons <- c(included_exons,paste(unique_including,sep=" ", collapse = " "))
        }else{
            included_exons <- c(included_exons,NA)
        }
        
    }
    
    
    ref <- as.data.frame(rbind(ref[c(1:4),],included_exons))
    print("may not work with minus stranded genes --> see comparison of refernces in adgrl2 in FIL")
    
    
    return(ref)
    
}

MakeMatchTable <- function(ref_comp_col, exons, locus){
    ref <- ref_comp_col
    
    ex_seq <- uncollapse(ref[5])[[1]]
    if(ex_seq[1]==""){
        ex_seq <- ex_seq[-1]
    }
    
    ex <- exons
    
    loc <- locus
    
    # get the sequence of the exons
    
    ex_seqs <- c()
    
    for(i in allrows(exons)){
        #print(i)
        # i <- 1
        coords <- c(exons[i,c(4)],exons[i,c(5)])
        
        coords <- coords-exons[1,4]+1
        
        
        loc_i <- s2c(loc$V1[2])
        ex_seqs <- c(ex_seqs,c2s(loc_i[coords[1]:coords[2]]))
        
        
    }
    
    ex <- cbind(ex,ex_seqs)
    x <- sapply(ex, is.factor)
    ex[x] <- lapply(ex[x], as.character)        
    
    
    match_table <- c()
    for(n in ex_seq){
        #print(n)
        #n <- ex_seq[1]
        n <- as.numeric(n)
        seq_n <- ex$ex_seqs[n]
        match_res <- matchPattern(pattern = seq_n, subject = ref_seq)
        match_table <- rbind(match_table, c(n,start(match_res),end(match_res)))
    }
    
    match_table <- as.data.frame(match_table)
    
    # check distances
    
    for(n in allrows(match_table)){
        # n <- 1
        if(n == 1){
            dist <- c(match_table$V2[n])
        }else{
            
            dist_n <- match_table$V2[n]-match_table$V3[(n-1)]
            dist <- c(dist,dist_n)
        }
        
    }
    match_table <- cbind(match_table,dist)
    
    colnames(match_table) <- c("exon","align_start","align_end","distance")
    
    
    name <- c2s(s2c(ref[1])[1:20])
    
    x <- c()
    for(n in allrows(match_table)){
        #print(n)
        x <- c(x,c(match_table[n,2]:match_table[n,3]))
    }
    
    #plot(x,main=name)
    
    plot_table <- c()
    for(n in 1:max(x)){
        if(n%in%x){
            plot_table <- rbind(plot_table,c(n,1))
        }else{
            plot_table <- rbind(plot_table,c(n,0))
        }
    }
    
    plot(plot_table[,2]~plot_table[,1], main = name)
    
    return(match_table)
    
}

RelativeTranscriptComparisonWithSD <- function(unique_transcripts = unique_transcripts, transcripts = transcripts, genename_length = 6){
    
    print("Careful, further pipeline with output of this functino NOT guaranteed, if you want to be safe use function without SD")
    
    #tr <- unique_transcripts
    tr <- transcripts
    tr_uniq <- unique_transcripts
    
    # get relative values
    tr_rel <- AddSamplesAndRelFPKMToTranscripts(tr_table = tr, genename_length = genename_length)
    
    # convers unique_transcripts
    freqs <- c()
    for(i in allrows(tr_uniq)){
        #i <- 1
        freqs <- c(freqs,length(uncollapse(tr_uniq$V2[i])[[1]]))
    }
    tr_uniq <- cbind(tr_uniq,freqs)
    tr_uniq <- tr_uniq[order(tr_uniq$freqs)[length(tr_uniq$freqs):1],-3]
    
    # now add coloumns to tr_uniq
    
    tr_uniq <- cbind(tr_uniq,NA,0,0)
    colnames(tr_uniq)[3] <- "all"
    colnames(tr_uniq)[4] <- "all_rel"
    colnames(tr_uniq)[5] <- "all_rel_sd"
    
    
    # how much coloumns do you need?
    
    tissues <- as.character(unique(tr_rel$type))
    
    len_new_cols <- length(tissues) * 3
    
    for(n in tissues){
        #n <- tissues[1]
        tr_uniq <- cbind(tr_uniq,NA,0,0)
        colnames(tr_uniq)[c((ncol(tr_uniq)-2):ncol(tr_uniq))] <- c(n,paste(n,"_rel",sep="",collapse=""),paste(n,"_rel_sd",sep="",collapse=""))
    }
    
    # fill it wiht values and tissue dependency
    
    for (i in allrows(tr_uniq)){
        # i <- 1
        samples <- uncollapse(tr_uniq$V2[i])[[1]]
        # add and in the end divide
        for(n in samples){
            #n <- samples[1]
            tr_tmp <- tr_rel[tr_rel$V1==n,]
            tissue_tmp <- tr_tmp$type[1]
            rel_tmp <- as.numeric(RemoveElements(char = uncollapse(tr_tmp$V7[1])[[1]][8],amount = 1))
            
            col_tmp <- which( colnames(tr_uniq)==tissue_tmp )
            tr_uniq[i,col_tmp] <- paste(tr_uniq[i,col_tmp],n)
            tr_uniq[i,(col_tmp+1)] <- tr_uniq[i,(col_tmp+1)]+rel_tmp
            
        }
        
        
    }
    
    indiv <- unique(tr_rel$affil)
    tr_uniq$all_rel <- as.numeric(tr_uniq$all_rel)
    
    for(i in allrows(tr_uniq)){
        # i <- 1
        tr_uniq$all_rel[i] <- sum(as.data.frame(tr_uniq[i,grepl(x = colnames(tr_uniq), pattern = "_rel")]))
        tr_uniq$all_rel[i] <- tr_uniq$all_rel[i]/length(indiv)
        
    }
    
    
    for(i in tissues){
        # i <- tissues[1]
        col_tmp <- which( colnames(tr_uniq)==i )+1
        
        indiv_tmp <- indiv[grepl(pattern = i, x = indiv)]
        replicates_tmp <- length(unique(indiv_tmp))
        
        tr_uniq[,col_tmp] <- tr_uniq[,col_tmp]/replicates_tmp
        
    }
    
    
    for(i in tissues){
        # i <- tissues[2]
        col_tmp <- which( colnames(tr_uniq)==i )
        for(n in allrows(tr_uniq)){
            # n <- 1
            if(is.na(tr_uniq[n,i])){
                # dont do anything
            }else{
                tr_uniq[n,i] <- RemoveElements(tr_uniq[n,i], amount = 3, start = "first")
            }
            
        }
    }
    
    tr_uniq_finished <- tr_uniq[,-3]
    
    # now calculate SD
    
    cols_sd <- length(tissues)+1
    cols_sd <- 1:cols_sd
    #cols_sd <- colnames(tr_uniq_finished)[cols_sd]
    cols_sd <- -1+(3*cols_sd)
    
    # how many samples are there?
    samples_all <- length(unique(x = tr_rel$affil))
    
    for(i in cols_sd){
        #i <- cols_sd[1]
        for(n in allrows(tr_uniq_finished)){
            #n <- 1
            list_tmp <- uncollapse(tr_uniq_finished[n,i])[[1]]
            values_tmp <- c()
            for(x in list_tmp){
                tr_tmp <- tr_rel[tr_rel$V1==x,]
                values_tmp <- c(values_tmp,as.numeric(RemoveElements(char = uncollapse(tr_tmp$V7[1])[[1]][8],amount = 1)))
            }
            while(length(values_tmp)<samples_all){
                values_tmp <- c(values_tmp,0)
            }
            tr_uniq_finished[n,i+2] <- sd(values_tmp)
        }
    }
    
    
    
    
    return(tr_uniq_finished)
    
}













