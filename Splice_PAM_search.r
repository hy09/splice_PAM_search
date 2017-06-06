#!/usr/bin/env Rscript

# Huayun Hou, April, 2017.
# This script takes a gtf file of a specific chromosome, and a file denoting the overlap between HGMD mutations and exon-intron junctions. For each mutation, it searches first the exon-intron junction region for PAM sequences which allows a cutting to eliminate the mutation ("mutation-eliminating cut"). In the cases where such a PAM sequence is found, it keeps searching into up to 1000bp into the intron for another PAM sequence which allows to restore 5'splice donor site ("Splice-site generating cut"). The script reports all the search result in an R object and outputs summary tables. 

# Load libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Biostrings))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(GenomicRanges))


# Read in args
args<-commandArgs(TRUE)
print(args)
chr <- args[1]
gtf_file <- args[2]
hgmd_ol <- args[3]

# Read in a list of PAM sequences to search for
PAMs <- scan("PAM_sequences.txt", what = "character")
# use the maximum length of the PAM sequence to determine search space for PAM sequences at junctions
maxPAMlength <- max(sapply(PAMs, nchar))
# the relative location of +1 site after cutting
junction_loc <- maxPAMlength + 4

# Define genome data. 
hg19_genome <- BSgenome.Hsapiens.UCSC.hg19

# Source functions used in this script.
source("exon_intron_PAM_search_functions.r")


# for testing
# reading in gtf file
#setwd("~/mdwilson/huayun/collabs/cohn/test")

#gtf_file <- "test_chr22.gtf"
#hgmd_ol <- "~/mdwilson/huayun/collabs/cohn/junctions/hg19_v3/Exon_intron_junction_HGMD_plus6_chr22.bed"
#gtf_file <- "lama2.gtf"
#gtf_file <- "~/mdwilson/huayun/collabs/cohn/gtf/hg19/gencode.v19.annotation.chr22.gtf"

# Read in the gtf file of a specific chromosome.
gtf <- read.table(gtf_file, fill = T, as.is=T, quote="\"", sep="\t")

# Read in a file of hgmd entries which are within 6bp from an intron-exon junction. 
hgmd <- read.table(hgmd_ol, fill=T, as.is=T, sep="\t", quote="")
# Reorder the hgmd entries table and convert it to GRanges.
hgmd <- hgmd[,c(14:20, 1:13, 21)]
names(hgmd)[c(1,2,3,4,6)] <- c("chr", "start", "end", "gene_id", "strand")
hgmd_gr <- makeGRangesFromDataFrame(hgmd, keep.extra.columns = TRUE)

# Extract only genes have hgmd overlap and entries annotated as "CDS" from the gtf file. 
gtf$gene_id <- get_txid(gtf$V9, type="gene_id")
gtf <- subset(gtf, V3 == "CDS" & gene_id %in% hgmd$gene_id)

# only look at coding sequence, requires GTF file to have the CDS as an annotated featured.

# Get a match between gene id and gene names 
name_match <- unique(data.frame(gene_id = get_txid(gtf$V9, "gene_id"),
                                gene_name = get_txid(gtf$V9, "gene_name")))
name_match_2 <- name_match$gene_name
names(name_match_2) <- name_match$gene_id
name_match <- name_match_2

# Split the gtf file into a list of data frames by transcript id. 
gtf$txid <- get_txid(gtf$V9)
gtf_list <- split(gtf, gtf$txid)
  
# Single exon transcripts are excluded
gtf_list <- gtf_list[sapply(gtf_list, function(x) nrow(x) > 1)]
  
# For each transcript, get their intronic regions, separately for "+" and "-" strand transcripts. 
intron_list <- lapply(gtf_list, function(x){
  if (x[1,7] == "+"){
    get_introns_plus(cds = x, PAMlength = maxPAMlength)
  } else {
    get_introns_neg(cds = x, PAMlength = maxPAMlength)
  }
})
  
# Collapse results to one data frame 
intron_df <- do.call("rbind", intron_list)
names(intron_df) <- c("chr", "start", "end", "strand", "gene_id", "transcript_id", "intron_length", "index")
# remove the large list
rm(intron_list)
  
# Collapse same introns from different transcripts of the same gene
intron_collapsed <- intron_df %>% 
  group_by(gene_id) %>%
  do(collapse_introns(.)) %>%
  mutate(length = as.integer(end) - as.integer(start) + 1) %>%
  filter(length > 2*(maxPAMlength + 3)) %>%
  as.data.frame()
# Remove the large data frame of all introns
rm(intron_df)

# Transfer intron_collapsed data frame into Granges
intron_collapsed$start <- as.numeric(intron_collapsed$start)
intron_collapsed$end <- as.numeric(intron_collapsed$end)
intron_collapsed_gr <- makeGRangesFromDataFrame(intron_collapsed, keep.extra.columns = TRUE)
  
# Intersect with hgmd data to get the introns that overlap with hgmd mutations
ol <- findOverlaps(query=hgmd_gr, subject=intron_collapsed_gr, type="within", select="all")
intron_collapsed_withMutations <- data.frame(intron_collapsed_gr[subjectHits(ol),], hgmd_gr[queryHits(ol),])
intron_collapsed_withMutations <- intron_collapsed_withMutations[,c(1:3,5:10,20,22,25,26,27,28,29,32)]
names(intron_collapsed_withMutations)[10:ncol(intron_collapsed_withMutations)] <- c("mut_pos", "HGMD_id","HGMD_type", "SNP_id", "gene_name", "gene_description", "mut_sequence", "relative_pos")


# For each intron sequence, first see if PAM sequences were detected at the right position near junctions. If a junction-cutting site is found, search up to 1000bp into the intronic regions for a second PAM sequence which allows a cutting that can restore 5' splice donor site. 
# currently only works for PAM sequence which allow cutting right at exon-intron junction (first nucleotide of introns).

# Print information 
print(paste("Analyzing", length(unique(hgmd$V4)), ":", length(unique(intron_collapsed_withMutations$HGMD_id)), "HGMD mutations on", chr))
# Remove duplicated entries
intron_collapsed_withMutations <- intron_collapsed_withMutations[sapply(1:nrow(intron_collapsed_withMutations), function(x) judge(intron_collapsed_withMutations[x,])),]
# sometimes a same mutation corresponds to different gene names, in those cases just remove duplicates.
intron_collapsed_withMutations <- intron_collapsed_withMutations[!duplicated(intron_collapsed_withMutations$HGMD_id),]
# Check if all the HGMD mutations are captured. 
if(nrow(intron_collapsed_withMutations) != length(unique(hgmd$V4))){
  stop("Numbers of unique HGMD ids do not match")
}

# Get all unique HGMD ids
HGMD_ids <- intron_collapsed_withMutations$HGMD_id
# Use the "get_seq" function to get results for PAM sequence searching in the exon-intron junction as well as intronic regions following the HGMD mutation. 
PAM_match_all_list <- sapply(HGMD_ids, function(x) get_seq(subset(intron_collapsed_withMutations, HGMD_id == x), hg19_genome))
names(PAM_match_all_list) <- HGMD_ids

# Remove the searches which result in no PAMs found to make a cutting at exon-intron junctions.
PAM_match_all_list_noNA <- PAM_match_all_list[!is.na(PAM_match_all_list)]
intron_count <- sapply(PAM_match_all_list_noNA, "[[", 3)
PAM_match_all_list_noNA_withintron <- PAM_match_all_list_noNA[names(intron_count[intron_count!=0])]
# A summary data frame of locations of each HGMD mutations, its corresponding gene name, relative position of the mutations to exon-intron junction and how many "mutation-eliminating cut (cut at exon-intron junctions)" and "splice-site generating cut" (cut within intron to restore a donor site) are found.
PAM_match_all_sum_df <- na.omit(do.call("rbind", lapply(PAM_match_all_list_noNA_withintron, output_sum)))

# Save all the search results in to an R object
output_rds <- paste0("Splice_PAMs_", chr, "_", Sys.Date(), ".rds")
saveRDS(PAM_match_all_list, file=output_rds)

# Save the summary table into a text file.
output_file <- paste0("Splice_PAMs_", chr, "_", Sys.Date(), ".txt")
write.table(PAM_match_all_sum_df, output_file, row.names=FALSE, col.names=TRUE, quote=F, sep="\t")

# make a summary


summary_df <- data.frame(HGMD_id = HGMD_ids,
                         HGMD_mut_pos = sapply(HGMD_ids, function(x) subset(intron_collapsed_withMutations, HGMD_id == x)$relative_pos) + 1,
                         junction_PAM = HGMD_ids %in% names(PAM_match_all_list_noNA),
                         intronic_PAM = HGMD_ids %in% names(PAM_match_all_list_noNA_withintron))

# Make a summary of number of HGMD mutations at each position (+1 means its at the first nucleotide of introns), the total number of "mutation-eliminating cut" and "splice-site generating cut" found for each case.
summary_df <- summary_df %>% 
  group_by(HGMD_mut_pos) %>%
  summarize(count = length(HGMD_mut_pos),
            junction_PAM = length(junction_PAM[junction_PAM]),
            intronic_PAM = length(intronic_PAM[intronic_PAM])) %>%
  as.data.frame()
write.table(summary_df, paste0("Splice_PAMs_", chr, "_", Sys.Date(), "_summary.txt"), row.names=FALSE, col.names=TRUE, quote=F, sep="\t")



