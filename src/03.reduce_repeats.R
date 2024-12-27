# This script aims to resolve overlaps in repeat annotations in a reduce manner if present
pacman::p_load(rtracklayer, Rsamtools, dplyr)

faFile = "~/storage/Data/cEle_wagos_eco/raw/genome/escherichia_coli.PRJNA526029.genomic.fa"
gffFile = "escherichia_coli.PRJNA526029.genomic.gff"
outFile = sub("genomic.gff", "disjoin_repeats.gff3", gffFile)
# import and format
gff = import(gffFile)
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(gff) = seqlevels(sInfo)
seqinfo(gff) = sInfo

# extract RNA classes
reps = gff[gff$type == "direct_repeat",]
tes = gff[grepl("transposase",gff$product),]
gff = c(reps, tes)
# Remove useless levels
gff$type = as.character(gff$type)
gff$type = ifelse(gff$type == "CDS", "repreat", gff$type)
gff$type = factor(gff$type)

# I expect a few overlap but in any case I will use reduce
gffList = split(gff, gff$type)

# remove overlaps if present
redList = list()
for (i in names(gffList)) {
  # extract by class
  tmp.cls = gffList[[i]]
  # wdt to chose just in case
  tmp.cls$wdt = width(tmp.cls)
  # Reduce 
  tmp.red = reduce(tmp.cls, with.revmap = T)
  # define uniques and duplicates
  uniques = tmp.red[unlist(lapply(tmp.red$revmap, length)) == 1,]
  dups = tmp.red[unlist(lapply(tmp.red$revmap, length)) > 1,]
  # For unique assign meteainfo directly
  mcols(uniques) = mcols(gff[c(uniques$revmap %>% unlist()),])
  # If duplicated (not trivial in bacteria)
  if (length(dups) >= 1) {
    # split each duplicated element
    tmp.lts = split(dups)
    rdup_lst = list()
    for (j in seq_along(tmp.lts)) {
      # Obtain the info of overlapping elements 
      tmp.dup = tmp.lts[[j]]
      tmp.hts = tmp.cls[c(unlist(tmp.dup$revmap)),]
      # Assign the meta info of the longest
      tmp.hts = tmp.hts[order(tmp.hts$wdt, decreasing = T),][1] %>% mcols()
      mcols(tmp.dup) = tmp.hts
      rdup_lst[[j]] = tmp.dup
    }
    # format, merge and save
    names(rdup_lst) = NULL 
    rdup_lst = do.call(c, rdup_lst)
    all_ranges =  c(uniques, rdup_lst)
    redList[[i]] = all_ranges
  } else{
    # else just save uniques
    redList[[i]] = uniques
  }
}
# Format and save
names(redList) = NULL
redRanges = do.call(c, redList)
export(redRanges, outFile, "gff3")



