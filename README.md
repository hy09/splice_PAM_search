# splice_PAM_search

These scripts were used for searching PAM sequences required for different Cas9 species. The aim is to use CRISPR to target splice-donor site mutations at +1 (first nucleotide in the intron after a exon-intron junction) to +6 positions. Specifically, two cut site are needed; one generates a cut which eliminates the splice-donor site mutation, the other one generates another cut within the following intronic region to restore the splice-donor site. 


### get_exon_intron_junction_gtfToBed.r

  - Given a gtf file of a specific chromosome, extract the positions of all exon-intron junctions.
  - The same junction sites for different transcripts of one gene are reported as one entry.
  
### Splice_PAM_search.r

  - Used for searching required PAM sequences.
  - Input files are tables of splice-donor site HGMD mutations and their overlapping genes; as well as gtf files.
 
### exon_intron_PAM_search_functions.r

  - Functions used by Splice_PAM_search.r script. 
