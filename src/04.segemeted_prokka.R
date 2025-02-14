setwd("~/MacStorage/CeMbio/ceMbio/tst")
# This script aims to produce a non-overlapping annotation in prokaryotic genomes 
# This script is based in Chow et al., 2019 segmented annotation form Cei
# however it was adapted for bacterial genomes with less TEs. 
# Also here I consider the gene structure in bacterial genomes
pacman::p_load(rtracklayer, Rsamtools, GenomicRanges, optparse)

# Command line options
option_list = list(
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Genome fasta file"),
  make_option(c("-p", "--gene_gff"), type = "character", default = NULL,
              help = "Disjoined gene annotation"),
  make_option(c("-r", "--repeats_gff"), type = "character", default = NULL,
              help = "Disjoined annotation of repeats"),
  make_option(c("-n", "--ncrna_gff"), type = "character", default = NULL,
              help = "Disjoined annotation of ncRNAs"),
  make_option(c("-o", "--out_gff"), type = "character", default = NULL,
              help = "Outfile name of the segmented annotation")
  
)

# Parse command line arguments
opt = parse_args(OptionParser(option_list = option_list))
# Check if input file is provided
# Validate input arguments
if (is.null(opt$fasta)) {
  stop("Error: No fasta file provided. Use --fasta to specify the input file.")
}
if (is.null(opt$gene_gff)) {
  stop("Error: No gene annotation provided. Use --gene_gff to specify the input file.")
}

if (is.null(opt$repeats_gff)) {
  stop("Error: No repeat annotation provided. Use --repeats_gff to specify the input file.")
}

if (is.null(opt$ncrna_gff)) {
  stop("Error: No ncRNA annotation provided. Use --ncrna_gff to specify the input file.")
}

if (is.null(opt$out_gff)) {
  stop("Error: No out file name provided. Use --out_gff to specify the output file.")
}

# Read files
genome = opt$fasta # assembly
proteinsFile = opt$gene_gff # from 01
repeatsFile = opt$repeats_gff # from 03
ncRNAFile = opt$ncrna_gff # from 02
outFile = opt$out_gff # outfile

# import files
proteins = import(proteinsFile)
repreats = import(repeatsFile)
ncrna = import(ncRNAFile)



# I need a list with all elements
allList = list(proteins, repreats, ncrna)
names(allList) = c("gene", "repeat", "ncrna")

# In bacteria genes are organized in operons, and most genes lack intron-exon structure
# Since this can cause some ranges issues, 01 script uses disjoint to assign overlapping bases
# However, after inspecting with IGV these annotations I noticed that relevant functions like 
# setdiff may fail, merging contiguous ranges. Bearing this in mind I decided to define a variable
# to artificially introduce a small space between genes. desired_shrink will help to shrink my 
# gene annotations by 2 bases in each edge. 
desired_shrink = 2

# The next section can easily out of a for but for control I'll it keep this way
strandedList = list()
for (i in names(allList)) {
  # first ncRNAs spited by type
  if (i == "ncrna") {
    tmp.class = allList[[i]]
    sInfo = seqinfo(scanFaIndex(genome))
    seqlevels(tmp.class) = seqlevels(sInfo)
    seqinfo(tmp.class) = sInfo
    ncList = split(tmp.class,tmp.class$type)
    nList = list()
    for (j in names(ncList)) {
      x = ncList[[j]]
      mcols(x) = NULL
      x$type = j
      x$class = paste0(x$type, "_S")
      x_as = x 
      strand(x_as) = ifelse(strand(x_as) == "+", "-", "+")
      x_as$class = sub("_S", "_As", x_as$class)
      n_nc = c(x, x_as)
      nList[[j]] = n_nc
    }
    names(nList) = NULL
    ncRanges =  do.call(c, nList)
    strandedList[[i]] = ncRanges
  }
  # Now genes
  if (i == "gene") {
    tmp.class = allList[[i]]
    sInfo = seqinfo(scanFaIndex(genome))
    seqlevels(tmp.class) = seqlevels(sInfo)
    seqinfo(tmp.class) = sInfo
    mcols(tmp.class) = NULL
    tmp.class$type = "CDS"
    tmp.class$class = "CDS_S"
    x_as = tmp.class
    strand(x_as) = ifelse(strand(x_as) == "+", "-", "+")
    x_as$class = sub("_S", "_As", x_as$class)
    prot_list = c(tmp.class, x_as)
    # Resize is a key element of this script since very contiguous ranges can
    # cause unexpected behaviors by functions like setdiff 
    prot_list = resize(prot_list, width = width(prot_list) - (2 * 2), fix = "center")
    strandedList[[i]] = prot_list
  }
  # Now repeats
  if (i == "repeat") {
    tmp.class = allList[[i]]
    sInfo = seqinfo(scanFaIndex(genome))
    seqlevels(tmp.class) = seqlevels(sInfo)
    seqinfo(tmp.class) = sInfo
    mcols(tmp.class) = NULL
    tmp.class$type = "repeat"
    tmp.class$class = "repeat_S"
    x_as = tmp.class
    strand(x_as) = ifelse(strand(x_as) == "+", "-", "+")
    x_as$class = sub("_S", "_As", x_as$class)
    rep_list = c(tmp.class, x_as)
    strandedList[[i]] = rep_list
  }
}
names(strandedList) = NULL 
allSegemets = do.call(c, strandedList)
allSegemetsByCalss = split(allSegemets, allSegemets$class)


ncFams =  names(allSegemetsByCalss)[!grepl("CDS|repeat",names(allSegemetsByCalss))]
# I'm going to define a hierarchical classification to solve overlaps
strandedCategories = c(ncFams, "repeat_S", "repeat_As","CDS_S", "CDS_As")

stopifnot(names(allSegemetsByCalss) %in% strandedCategories)

strandedList = allSegemetsByCalss[strandedCategories]
allWidth = do.call("sum",sapply(c(strandedList), width))

# Start with one, cycle through the others
print("Removing overlaps, according to sorted types")
allSegments <- strandedList[[1]]
allSegments$class <- names(strandedList)[1]
for (type in names(strandedList)[-1]) {
  print(type)
  allSegments[order(seqnames(allSegments), start(allSegments)),]
  trimmedSegs <- setdiff(strandedList[[type]], allSegments)
  if (length(trimmedSegs) > 0) {
    trimmedSegs$class <- type
    allSegments <- c(allSegments, trimmedSegs)
  }
}
allSegments$type = NULL

# Calculate difference in sum of all widths after overlap reduction
remWidth <- round((allWidth-sum(width(allSegments)))/allWidth*100,2)
print(paste0("Removed ",remWidth,"% of annotated bases due to overlaps in sorted types"))


# Calculate all intergenic spaces
allGaps <- gaps(reduce(allSegments, ignore.strand=TRUE))
allGaps <- allGaps[strand(allGaps) == "*"]
allGaps$class <- "intergenic"

# Calculate the ranges that were removed due to any overlap
allSegments = c(allSegments, allGaps)
# Check genome sizes
allTypes <- sapply(split(width(allSegments), allSegments$class), sum)

# multiply *2 those that are on the '*' strand
for (type in c("intergenic")) {
  allTypes[type] <- allTypes[type] * 2
}
allTypes = allTypes[!is.na(allTypes)]

# Check point for genome coverage
print(paste("Annotated genome size=",sum(allTypes)/1e6/2,"Mb"))
print(paste("Real genome size=",sum(seqlengths(allGaps))/1e6,"Mb"))
ann_bases = allTypes/1e6/2
# Format and save
allSegments = allSegments[order(seqnames(allSegments), start(allSegments)),]
allSegments$ID = paste0(allSegments$class, ":", 1:length(allSegments))
export(allSegments, outFile, "gff3")











