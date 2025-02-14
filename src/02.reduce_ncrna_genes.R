# This script aims to resolve overlaps in ncRNA annotations in a reduce manner if present
pacman::p_load(rtracklayer, Rsamtools, dplyr,optparse)

# Command line options
option_list = list(
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Genome fasta file"),
  make_option(c("-g", "--gff"), type = "character", default = NULL,
              help = "Full genome annotation (typically from prokka)"),
  make_option(c("-o", "--out_gff"), type = "character", default = NULL,
              help = "Out gff file name")
)

# Parse command line arguments
opt = parse_args(OptionParser(option_list = option_list))
# Check if input file is provided
# Validate input arguments
if (is.null(opt$fasta)) {
  stop("Error: No fasta file provided. Use --fasta to specify the input file.")
}
if (is.null(opt$gff)) {
  stop("Error: No annotation provided. Use --gff to specify the input file.")
}
if (is.null(opt$out_gff)) {
  stop("Error: No out file name provided. Use --out_gff to specify the output file.")
}

# Pass input variables
faFile = opt$fasta
gffFile = opt$gff
outFile = opt$out_gff

# import and format
gff = import(gffFile)
sInfo = seqinfo(scanFaIndex(faFile))
seqlevels(gff) = seqlevels(sInfo)
seqinfo(gff) = sInfo

# extract RNA classes
gff = gff[gff$type %in% c("rRNA","tRNA", "ncRNA", "riboswitch","antisense_RNA", 
                          "SRP_RNA", "tmRNA", "sequence_feature"),]


# Remove useless levels
gff$type = as.character(gff$type)
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
