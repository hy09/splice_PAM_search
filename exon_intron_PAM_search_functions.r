# Functions used by "Splice_PAM_search.r" script.

get_txid <- function(anno, type="transcript_id"){
  # Function to get transcript id from the annotation field.
    anno_list <- unlist(strsplit(anno, "; "))
    return(gsub(paste0(type, " "), "", anno_list[grepl(type, anno_list)]))
}
# }

# Function to collapse same introns from different transcripts of the same gene
# Preserve the smallest intron size for collapsed introns
collapse_introns <- function(gene_df){
  # group by chr and start, choose the smallest end (for plus strand genes) to resolve cases for overlapping exons
 if(gene_df[1,4] == "+"){
    gene_df <- gene_df %>% 
      group_by(chr, start) %>%
      summarize(end = min(end), strand = unique(strand), gene_id = unique(gene_id), transcript_id = paste(sort(transcript_id), collapse=";"), intron_length = min(intron_length)) %>%
      as.data.frame()
    gene_df$index <- paste0("intron_", 1:nrow(gene_df))
  } else{
    gene_df <- gene_df %>% 
      group_by(chr, end) %>%
      summarize(start = max(start), strand = unique(strand), gene_id = unique(gene_id), transcript_id = paste(sort(transcript_id), collapse=";"), intron_length = min(intron_length)) %>%
      as.data.frame()
    gene_df$index <- rev(paste0("intron_", 1:nrow(gene_df)))
  }
  return(gene_df)
}


get_seq <- function(line, genome=hg19_genome, consensus_length = 2){
  # Main function for PAM sequence search. Parameters:
  # line: each line denoting intronic regions with an HGMD mutation as input.
  # consensus_length: how many bps are required to restore 5' splice donor site. Currently only considering the first 2 nucleotides, usually a "GT". 
  
  # Get the sequence of exon-extend + intron-extend, on the strand the transcript is on.
  intron_length <- as.numeric(line[1,"intron_length"])
  # The position of mutation, relative to exon-intron junction
  relative_mut_pos <- as.numeric(line[1, "relative_pos"]) + 1
  # Get the sequence 
  seq <- getSeq(genome, names = as.character(line[1,"seqnames"]), start = as.integer(line[1,"start"]), end = as.integer(line[1,"end"]), strand = as.character(line[1,"strand"]))

  # search the 2 * (maxPAMlength + 3) + extra bases for +N sites region around exon-intron junction first to look for all instances of PAM sequences
   junction_region <- seq[1:(2*(maxPAMlength + 3) + relative_mut_pos - 1)]
  
  # search for each PAM at junctions 
  PAMmatch_list <- lapply(PAMs, function(x) PAM_match(junction_region, x, intron_length, relative_mut_pos = relative_mut_pos))
  PAMmatch_junction_df <- unique(na.omit(do.call("rbind", PAMmatch_list)))
  nrow(PAMmatch_junction_df)
  
  # if junction PAM is found, keep searching
  # TO DO: deal with cutting at +1 or +2 sites
  # get the 2 nucleotides at the exon-intronic junction
  # TO DO: make this a flexible length
  junction_plus1 <- maxPAMlength + 4
  junction_nuc_plus1 <- as.character(seq[junction_plus1:(junction_plus1 + consensus_length -1)])
  
  if (nrow(PAMmatch_junction_df) > 0){
    # if any junction cuts are found, keep looking into the introns. Find all the PAMs within the intron.
    PAMmatch_all_df <- unique(na.omit(do.call("rbind", lapply(PAMs, function(x) PAM_match(seq, x, intron_length, all = TRUE, relative_mut_pos = relative_mut_pos)))))
    # Remove the PAMs results in cutting at the very end of the introns
    # TO DO: allow flexible searching in these cases.
    PAMmatch_all_df_fil <- subset(PAMmatch_all_df, cutsite < (nchar(seq) - consensus_length)) 
    # Get a list of where are junction cuts found
    junction_cut_pos <- unique(PAMmatch_junction_df$cutsite) - (maxPAMlength + 3)
    if(length(setdiff(junction_cut_pos, c(1, 2))) == 0){
      # if the cut is at +1 or +2, search if the intronic cut could restore the "GT".
      # TO DO : report the sequence after the cut
      # Get the nucleotides follwoing the cut sites in introns.
      cutsite_nucs <- sapply(PAMmatch_all_df_fil$cutsite, function(x) as.character(seq[x:(x + consensus_length - 1)]))
      # initiate logical vectors to store judgement for +1 or +2 position nucleotide matches
      PAM_true_1 <- PAM_true_2 <- rep(FALSE, length(cutsite_nucs))
      if(1 %in% junction_cut_pos){
        # if some of the junction cuts are made at +1 position, need to make sure both +1 and +2 nucleotides are restored.
        PAM_true_1 <- cutsite_nucs == junction_nuc_plus1
      } else if(2 %in% junction_cut_pos){
        # if the cutting is at the +2 position, only need to make sure the +2 nucleotide is restored.
        PAM_true_2 <- unname(sapply(cutsite_nucs, function(x) substr(x, 1, 1))) == substr(junction_nuc_plus1,1,1)
      }
      # Get a list of PAMs whose guided cutting will restore the "GT" at exon-intron junction.
       PAM_true <- PAM_true_1 | PAM_true_2 
       PAMmatch_all_df_fil <- PAMmatch_all_df_fil[PAM_true,]
    } 
        
    if (nrow(PAMmatch_all_df_fil) > 0){
      # If at least one PAM is found to allow a "splice-site generating cut" within intron, report the result. For now, only return junction hits and number of intronic cuts as the specific information may be too extensive and can be easily re-generated if focusing on only on intron.

      output_PAMs <- PAMmatch_to_output(line, PAMmatch_junction_df)
      nIntronPAMs <- nrow(PAMmatch_all_df_fil)
      uniqIntronPAMs <- length(unique(PAMmatch_all_df_fil$cutsite))
      output <- list(junctionPAMs = output_PAMs, 
                          intron_info = line,
                          nIntronPAMs = nIntronPAMs,
                          uniqueIntronPAMs = uniqIntronPAMs)
    } else {
      # If no PAMs were found within the intron. Report 0 for "nIntronPAMs" but still report the PAMs found at junctions. 
      output <- list(junctionPAMs = PAMmatch_junction_df,
                     intron_info = line,
                     nIntronPAMs = 0,
                     uniqueIntronPAMs = 0)
    }
  } else{
    # If no PAMs were found to make "mutation-eliminating cut" at junctions, report NA.
    output <- NA
  }
  return(output)
}

PAMmatch_to_output <- function(intron_line, potential_PAMs){
  # this function takes the intron_line and a data frame of PAM hits and generates the output dataframe.
    if(as.character(intron_line["strand"])=="+"){
      potential_PAMs$genomic_start <- sapply(potential_PAMs$start, function(x) to_genomic(intron_line, x))
      potential_PAMs$genomic_end <- sapply(potential_PAMs$end, function(x) to_genomic(intron_line, x))
    } else {
      potential_PAMs$genomic_end <- sapply(potential_PAMs$start, function(x) to_genomic(intron_line, x))
      potential_PAMs$genomic_start <- sapply(potential_PAMs$end, function(x) to_genomic(intron_line, x))
    }
    
    potential_PAMs$cutsite <- potential_PAMs$cutsite - (maxPAMlength + 3)
    potential_PAMs$genomic_strand <- unname(sapply(as.character(potential_PAMs$strand), function(x) ifelse(x=="opp", reverse_strand(as.character(intron_line[1,"strand"])), as.character(intron_line[1,"strand"]))))
    return(potential_PAMs)
}



reverse_strand <- function(x){
  ifelse(x=="+", "-", "+")
}

to_genomic <- function(intron_line, relative_pos){
  # Function to translate relative position on intron sequence to genomic position
  g_start <- as.integer(as.character(intron_line["start"]))
  g_end <- as.integer(as.character(intron_line["end"]))
  g_strand <- as.character(intron_line["strand"])
  if(g_strand == "+"){
    g_pos <- g_start + relative_pos -1
  } else{
    g_pos <- g_end - relative_pos + 1 
  }
}

get_cut_site <- function(PAMmatch){
  # Function to get the cut site for a PAM sequence match use the first genomic location after the cut as index
  if(PAMmatch[3] == "opp"){
    cut_site <- as.integer(as.character(PAMmatch[2])) + 4
  } else {
    cut_site <- as.integer(as.character(PAMmatch[1])) - 3
  }
  return(cut_site)
}


get_gDNA <- function(PAMmatch, gDNAlength=20, seq){
  # Function to get the relative position of gDNA on the given sequence 
  if(PAMmatch[3] == "opp"){
    gDNAstart = as.integer(PAMmatch[2]) + 1
    gDNAend = as.integer(PAMmatch[2]) + gDNAlength
    gDNA = reverseComplement(seq[gDNAstart:gDNAend])
  } else{
    gDNAstart = as.integer(PAMmatch[1]) - gDNAlength
    gDNAend = as.integer(PAMmatch[1]) - 1
    gDNA = seq[gDNAstart:gDNAend]
  }
  return(as.character(gDNA))
}



get_introns_plus <- function(cds, PAMlength = 8, intron_extend = 1000){
  # function to get introns (intervals between exons), works for plus strand
  # 1) Extend into exonic regions by certain bps, allowing cut to occur at +1 region.
  # 2) Extend into intronic regions, making sure a. fixed maximum length; b. cut site wont occur within 50bp before next exon.
  
  # order data frame of exons by genomic positions
  cds <- cds[with(cds, order(V4, V5)),]
  
  # initiate data frame for intron
  introns <- as.data.frame(matrix(nrow=nrow(cds) -1,ncol=7))
  # loop through all exons to get intron intervals
  for (i in 1:(nrow(cds) - 1)){
    line1 <- cds[i,] # last exon
    line2 <- cds[i+1,] # next exon 
    intron_start_pos <- line1[1,5] + 1 # position of the first nucleotide of intron 
    intron_end_pos <- line2[1,4] - 1 # position of the last nucleotide of intron, check cut site position later
    # position of the last nucleotide of intron, minus 50bp to not interfer with donar site regulation 
    intron_length <- intron_end_pos - intron_start_pos + 1 # get intron length
    intron_end_pos <- intron_start_pos + min(intron_extend, intron_length) - 1 # if intron is longer than expected search space, set new intron end. 
    intron_start_extend_pos <- intron_start_pos - (PAMlength + 3) # this is the start site for PAM sequence searching in exon
    chr <- line1[1,1]
    strand <- line1[1,7]
    newline <- c(chr, intron_start_extend_pos, intron_end_pos, strand, get_txid(line1[1,9], "gene_id"), get_txid(line1[1,9], 
    introns[i,] <- newline
  }
  introns$index <- paste0("intron", 1:nrow(introns)) # index exons
  return(introns)
}


get_introns_neg <- function(cds, PAMlength = 8, intron_extend = 1000){
  # function to get introns (intervals between exons), works for negative strand
  # order data frame of exons by genomic positions. Note that for negative strand genes, exons order will be revsersed.
 
  cds <- cds[with(cds, order(V4, V5)),]

  # initiate data frame for intron
  introns <- as.data.frame(matrix(nrow=nrow(cds) -1,ncol=7))
  # loop through all exons to get intron intervals
  for (i in 1:(nrow(cds) - 1)){
    line1 <- cds[i,] # last exon
    line2 <- cds[i+1,] # next exon 
    intron_end_pos <- line1[1,5] + 1 # position of the end of intron
    intron_start_pos <- line2[1,4] - 1 # position of the last nucleotide of intron 
    intron_length <- intron_start_pos - intron_end_pos # get intron length
    intron_end_pos <- intron_start_pos - (min(intron_extend, intron_length) -1) # if intron is longer than expected search space, set new intron end. 
    intron_start_extend_pos <- intron_start_pos + (PAMlength + 3) # this is the start site for PAM sequence searching in exon
    chr <- line1[1,1]
    strand <- line1[1,7]
    newline <- c(chr, intron_end_pos, intron_start_extend_pos, strand, get_txid(line1[1,9], "gene_id"), get_txid(line1[1,9], "transcript_id"), intron_length)
    introns[i,] <- newline
  }
  introns$index <- rev(paste0("intron_", 1:nrow(introns))) # index intron
  return(introns)
}


PAM_match <- function(seq, PAM = "NNGRRT", intron_length, all = FALSE, relative_mut_pos){
  # This is the basic function searching one sequence for one PAM sequence. And if match, report matched sequence, cut site (relative to the sequence given) 
  seq <- DNAString(seq)
  # find PAM sequences match on the same strand as the gene
  same_hits <- matchPattern(PAM, seq, fixed = FALSE)
  # find PAM sequences macth on the opposite strand 
  opp_hits <- matchPattern(reverseComplement(DNAString(PAM)), seq, fixed = FALSE)
  PAMlength <- nchar(PAM)
  # process hits on the same strand
  if(length(same_hits) != 0) {
    same_df <- data.frame(start=start(same_hits), end=end(same_hits), strand="same", matched=as.character(same_hits))
  } else {
    same_df <- data.frame()
  }
  if(length(opp_hits) != 0) {
    opp_df <- data.frame(start=start(opp_hits), end=end(opp_hits), strand="opp", matched=as.character(reverseComplement(opp_hits)))
  } else {
    opp_df <- data.frame()
  }
  if(nrow(same_df) == 0 & nrow(opp_df) == 0){
    # If both strands had no hits, just return the junction search results. As there will be no need to search the whole sequence
      return(NA)
  } else{
    # If there are hits, assemble a data frame of all hits. 
    hit_df <- rbind(same_df, opp_df)
    # Get cut site for each hit_df
    hit_df$cutsite <- apply(hit_df, 1, get_cut_site)
    # Make sure cuts do not happen within 50bp of next exon
    hit_df <- subset(hit_df, cutsite <= (intron_length - 50))
    # If any hits can result in cut site at junction (+1)
    # TO DO: modify to include possibility of cutting at +2 
    junction_hit_df <- subset(hit_df, cutsite %in% c(junction_loc:(junction_loc + relative_mut_pos - 1)))
  
    # If any PAMs were found near junction, collapse the matched patterns
  
  if (!all){
    # If the function call was to check only PAM match at junctions
    if(nrow(junction_hit_df) > 0){
      return(junction_hit_df)
    } else {
      return(NA)
    }
  } else{
    # Only report cuts that are within introns and are out junction
    return(subset(hit_df, cutsite > (maxPAMlength + 3 + relative_mut_pos) & cutsite < (nchar(seq) -2)))
  }
}
  
# output_sum <- function(output_list){
#   HGMD_info <- output_list$intron_info[,c("seqnames", "mut_pos", "HGMD_id", "HGMD_type","gene_name","gene_description", "relative_pos")]
#   HGMD_info$relative_pos <- HGMD_info$relative_pos + 1
#   cut_sum <- table(output_list$junctionPAMs$cutsite)
#   HGMD_info$junctionPAMs <- paste(names(cut_sum), unname(cut_sum), sep=":", collapse=",")
#   HGMD_info$nIntronPAMs <- output_list$nIntronPAMs
#   return(HGMD_info)
# }
  
