# This script aims to remove overlapping bases based in prokaryotic gene annotations
pacman::p_load(rtracklayer, Rsamtools, dplyr)

# As input files I need the fasta and gff 
# For this script I will use E. coli OP50
# gff from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/355/015/GCF_004355015.1_ASM435501v1/

faFile = "~/storage/Data/cEle_wagos_eco/raw/genome/escherichia_coli.PRJNA526029.genomic.fa"
gffFile = "escherichia_coli.PRJNA526029.genomic.gff"
outFile = sub("genomic.gff", "disjon_annotation.gff3", gffFile)

# import and format
gff = import(gffFile)
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(gff) = seqlevels(sInfo)
seqinfo(gff) = sInfo

# I'm going to remove redundant annotations
# This is only for the E. coli genome
gff = gff[!grepl("gene|pseudogene|exon|region",gff$type),]
gff = gff[gff$type %in%  c("five_prime_UTR", "CDS", "three_prime_UTR", "three_prime_utr", "five_prime_utr"),]
# I may need to amend this for other gene annotations
gff = gff[!grepl("transposase", gff$product),]

# Now I would like to split to resolve 

# Find overlaps
overlaps = findOverlaps(gff, gff)
hits = overlaps[queryHits(overlaps) != subjectHits(overlaps),]
hits = hits[queryHits(hits) < subjectHits(hits),]

# Resolve overlaps
# Add lengths
gff$length = width(gff)

# Resolve overlaps
for (i in seq_along(hits)) {
  # Get overlapping pair
  overlap_pair = hits[i]
  # Extract the ranges
  overlapping_ranges = gff[c(queryHits(overlap_pair), subjectHits(overlap_pair))]
  # Ensure they are on the same strand
  if (length(unique(strand(overlapping_ranges))) <= 1) {
    # Order by length
    overlapping_ranges$length = width(overlapping_ranges)
    overlapping_ranges = overlapping_ranges[order(overlapping_ranges$length, decreasing = TRUE)]
    # Disjoin to split into unique segments
    disjointed = disjoin(overlapping_ranges, with.revmap = TRUE)
    # extact the shotest
    disjointed$poss = disjointed$revmap %>% paste(collapse = ",")
    short = disjointed[disjointed$poss == "2",]
    short$poss = NULL
    # resize  
    resized = GenomicRanges::setdiff(disjointed, short)
    resized$revmap = disjointed[disjointed$poss == "1",]$revmap
    resized = c(short, resized)
    # reassigned meta info based on "coverage"
    rs_list = list()
    for (j in seq_along(resized)) {
      tmp.gr = resized[j,]
      sub.tmp = subsetByOverlaps(gff, tmp.gr)
      sub_list =  split(sub.tmp)
      # Count coverage between overlaping regions
      sub.tmp$inter_length = lapply(sub_list, function(x){GenomicRanges::intersect(x, tmp.gr) %>% width()}) %>% unlist()
      # Order based on amount of covered regions 
      sub.tmp = sub.tmp[order(sub.tmp$inter_length, decreasing = T),]
      longer_gr = sub.tmp[1,]
      # Reasign based on the longest one
      mcols(tmp.gr) = mcols(longer_gr)
      rs_list[[j]] = tmp.gr
    }
    names(rs_list) = NULL
    resized = do.call(c, rs_list)
    resized$inter_length = NULL
    # Update the original GRanges object
    gff[c(queryHits(overlap_pair), subjectHits(overlap_pair))] = resized
  } else{
    # I will allow overlaps present on diffrent strand
    next()
  }
}
# Remove useless cols
gff$length = NULL
# And export
export(gff, outFile, "gff3")



