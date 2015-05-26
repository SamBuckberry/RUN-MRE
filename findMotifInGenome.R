## inspired by http://www.niravmalani.org/finding-motif-of-interest-in-a-genome-using-r/

library(Biostrings)
library(ggbio)
library(GenomicRanges)

## define the motif of interest ##
motif <- DNAStringSet(c("GGCC", "CCGG"))

## get the target sequence/genome of interest ##
target <- readDNAStringSet("~/genomes/Rattus_norvegicus/UCSC/rn5/Sequence/WholeGenomeFasta/genome.fa")
chroms <- names(target)

target <- target[[1]]  # extract the first sequence to make it a DNAString object

#### matchPattern(): long & roundabout ####
hits1 <- matchPattern(motif[[1]], target, fixed = TRUE)
hits2 <- matchPattern(motif[[2]], target, fixed = TRUE)

#### matchPattern(): short & simple ####
hits.all <- sapply(motif, matchPattern, subject = target, fixed = TRUE)

#### test to see if two methods shown above produce identical results ####
identical(hits.all[[1]], hits1)

sapply(hits.all, head)  # show the top 6 hits for all motifs

sapply(hits.all, length)  # show number of hits per motifs

#### make GRanges object out of list of all motif hits ####

## easy & straight forward way ##
hits.gr <- GRanges()  # initiate an empty GRanges object
for (x in hits.all) {
        hits.gr <- c(hits.gr, GRanges(seqnames = chroms[1], ranges(x)))  # append new GRanges to old
}
hits.gr


## add motif sequences to the GRanges object to differentiate hits per row
hits.gr$motif <- rep(as.character(motif), sapply(hits.all, length))

## add found motif sequences to the GRanges object
hits.gr$motif.found <- unlist(lapply(hits.all, as.character))
hits.gr

## add length of chromosome to GRanges object
chromLength <- length(target)
names(chromLength) <- as.character(chroms[1])
seqlengths(hits.gr) <- chromLength

# Add strand info
strand(hits.gr) <- "+"

## how many different motif variations found? ## Only useful if ambiguous char used
table(hits.gr$motif.found)

# Visualise the results

bamPath <- "~/tina_data/Project_JBCR_External/Sample_C583LANXX-7-Rat_GBS/allRatAlign.BAM.bam"
indexPath <- "~/tina_data/Project_JBCR_External/Sample_C583LANXX-7-Rat_GBS/allRatAlign.BAM.bam"

# Count the number of overlaps in the bam file
bam <- BamFile(file="~/tina_data/Project_JBCR_External/Sample_C583LANXX-7-Rat_GBS/allRatAlign.BAM.bam",
               asMates=TRUE, index="", yieldSize=1e5)

bamRecords <- scanBam(file=bamPath, index=bamPath, isProperPair=TRUE)

overlaps <- summarizeOverlaps(features=GRangesList(hits.gr), reads=bam, mode="union", singleEnd=FALSE, ignore.strand=TRUE)


