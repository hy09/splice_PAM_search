#!/usr/bin/env Rscript
# Created by Huayun Hou, April, 2017.  
# This script takes a gtf file of a specific chromosome and extracts all annotated exon-intron junctions.
# USAGE : Rscript get_exon_intron_junction_gtfToBed.r annotation.chr1.gtf chr1
# "annotation.chr1.gtf" referes to a gtf file of a specific chromosome. "chr1" denotes the chromosome name. 

# Load libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

# Read in args
args<-commandArgs(TRUE)
print(args)
gtf_file <- args[1]
chr <- args[2]

# Disable scientific notations
options(scipen = 999)

# Read in the gtf file 
gtf <- read.table(gtf_file, fill = T, as.is=T, quote="\"", sep="\t")

# Extract only features labeled as "CDS"
gtf <- subset(gtf, V3 == "CDS")

get_txid <- function(anno, type="transcript_id"){
  # Function to get transcript id from the annotation field
  anno_list <- unlist(strsplit(anno, "; "))
  return(gsub(paste0(type, " "), "", anno_list[grepl(type, anno_list)]))
}


get_introns <- function(cds){
  # Function to get introns (defined as intervals between exons of one transcript).
  # Input is a data frame of all "CDS" of a transcript.
  # Order data frame of exons by genomic positions. Note that for negative strand genes, exons order will be revsersed.
  cds <- cds[with(cds, order(V4, V5)),]
  
  # Initiate a data frame for intron
  introns <- as.data.frame(matrix(nrow=nrow(cds) -1,ncol=6))
  # Loop through all exons to get intron intervals
  for (i in 1:(nrow(cds) - 1)){
    line1 <- cds[i,] # last exon
    line2 <- cds[i+1,] # next exon 
    #intron_end_pos <- line1[1,5] + 1
    strand <- as.character(line1[1,7])
    txid <- as.character(line1[1, "txid"])
    gene_id <- as.character(line1[1, "gene_id"])
    if(as.character(strand) == "+"){
      intron_start_pos <- line1[1,5] + 1
    } else{
      intron_start_pos <- line2[1,4] - 1 # position of the last nucleotide of intron 
    }
    
    # Add gene id and transcript id information. 
    newline <- c(line1[1,1], intron_start_pos, intron_start_pos + 1, strand, gene_id, txid)
    # Write the new line into the intron data frame. 
    introns[i,] <- newline
  }
  if(as.character(strand) == "+"){
    introns$index <- paste0("intron_", 1:nrow(introns)) # index intron
  } else{
    introns$index <- rev(paste0("intron_", 1:nrow(introns))) # index intron, if negative strand, reverse the indices. 
  }
  return(introns)
}


collapse_gene <- function(gene_table){
  # Function to collapse exon-intron junctions from multiple transcripts for one gene. The transcript information is kept. 
  gene_collapsed <- gene_table %>%
    group_by(chr,start, end, strand, gene_id) %>%
    summarize(txid = paste(unique(txid), collapse=";"), 
              index = paste(unique(index), collapse=";")) %>%
    as.data.frame()
  return(gene_collapsed)
}

write_junctions <- function(gtf = gtf, chr="chr1"){
  # Function to write out a bed file of extracted intron-exon junctions. 
  gtf$txid <- get_txid(gtf$V9)
  gtf$gene_id <- get_txid(gtf$V9, type = "gene_id")
  gtf$tx_status <- get_txid(gtf$V9, type = "transcript_status")
  # filter out putative transcripts. Only keeping annotations annotated as "KNOWN".
  gtf <- subset(gtf, tx_status == "KNOWN")
  gtf_list <- split(gtf, gtf$txid)
  gtf_list <- gtf_list[sapply(gtf_list, function(x) nrow(x) > 1)]
  
  # For each transcript, get their intronic regions
  intron_list <- lapply(gtf_list, get_introns)
  
  # Merge results to one data frame 
  intron_df <- do.call("rbind", intron_list)
  names(intron_df) <- c("chr", "start", "end", "strand", "gene_id", "txid", "index")
  
  # For each gene, collapse the same intron-exon junctions from multiple transcripts. 
  gene_df <- do.call("rbind", lapply(split(intron_df, intron_df$gene_id), collapse_gene))
  # Reorder the data frame. 
  gene_df <- gene_df[, c("chr", "start", "end", "gene_id", "txid", "strand", "index")]
  gene_df <- gene_df[with(gene_df, order(start, end)),]
  # Create output file name
  output_file <- paste0("Exon_intron_junctions_", chr, ".bed")
  # Out put the data frame. 
  write.table(gene_df, output_file, row.names=FALSE, col.names=FALSE, quote=F, sep="\t")
}

write_junctions(gtf, chr=chr)

