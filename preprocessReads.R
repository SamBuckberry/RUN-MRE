library(ShortRead)
setwd(dir="~/tina_data/Project_JBCR_External/Sample_C583LANXX-7-Rat_GBS/")

# Specify fastq files
fls <- list.files(path="./",pattern="*_001.fastq.gz$", full.names=TRUE)

#### Load the barcodes file
barcodes <- read.table(file="barcodes.txt")
barcodes$motif <- "CGG"
barcodes$fullMotif <- paste(barcodes$V1, barcodes$motif, sep="")
barcodes$length <- nchar(barcodes$fullMotif)

##### Create the barcode filter ######

# Search by streaming fastq data
barcodeFinder <- function(fastqPath, barcode){
        
        # Setup the fastq streamer
        f <- FastqStreamer(con=fastqPath, n=1e7, verbose=TRUE)
        
        # Set function to apply filter to stream of reads
        while (length(fq <- yield(f))) {
                ## do work on reads here ###
                # setup the fastq accessor
                fq <- sread(fq)
                
                # Get the beginning of the reads to search for barcode
                fqSub <- substring(fq, first=1, last=nchar(barcode)+3)
                # Search for the barcode in the substring of the read and count
                count <- vcountPattern(pattern=barcode, subject=fqSub, fixed=TRUE,
                                       max.mismatch=0, with.indels=FALSE)
                # transform counts to logical
                match <- ifelse(test=count == 0, yes=FALSE, no=TRUE)
                
        }
        # Close the fastq streamer connection
        close(f)
        
        # Return the logical vector of barcode matches
        return(match)
}

# Run the filter over all barcodes and create a logical data.frame
t <- sapply(X=barcodes$fullMotif, FUN=barcodeFinder, fastqPath=fls[1])

# Count the number of different barcodes in each read.
tSums <- rowSums(t)
# Look at which reads contain more than one barcode
table(tSums)


#### Remove reads with ambiguous or no barcode in 5' leading sequence-----------
# Logical vector of reads with no barcode found
noBarcode <- tSums == 0
# Logical vector of more that one barcode
ambiguous <- tSums > 1
# Logical vector of clear barcode reads
goodReads <- tSums == 1

# Setup the simple function for filterFastq for noBarcode OR ambiguous reads
filterBad <- function(x){
        x[noBarcode | ambiguous]
}

# Extract bad reads output to new fastq file
filterFastq(files=c(fls[1], fls[2]),
            destinations=c("bad.R1.fastq.gz", "bad.R2.fastq.gz"),
            compress=TRUE, filter=filterBad, yieldSize=1e7)

# Proportion of bad reads
propBad <- sum(noBarcode, ambiguous)/nrow(t)
message(paste("Proprotion of bad reads = ", propBad, sep=""))

# Extract the good reads to a new file
filterGood <- function(x){
        x[goodReads]
}

filterFastq(files=c(fls[1], fls[2]), 
            destinations=c("good.R1.fastq.gz", "good.R2.fastq.gz"),
            filter=filterGood, # Only keep the reads with barcodes
            compress=TRUE, yieldSize=1e7)

#### Sepeate into seperate files based on barcode ----------------------

#barcode <- as.character(barcodes$V1[1])
#fullMotif <- barcodes$fullMotif[1]
R1 <- "good.R1.fastq.gz"
R2 <- "good.R2.fastq.gz"

demultiplexSamples <- function(R1, R2, barcode, fullMotif){
        
        # Setup the destination files
        outR1 <- paste(barcode, ".R1.fastq.gz", sep="")
        outR2 <- paste(barcode, ".R2.fastq.gz", sep="")
        
        # Determine which reads feature the barcode of interest
        readMatch <- barcodeFinder(fastqPath=R1, barcode=fullMotif)
        
        # Extract the good reads to a new file
        filterMatch <- function(x){
                x[readMatch]
        }
        
        # Write the reads to a new file with barcode as filename
        filterFastq(files=c(R1, R2), 
                    destinations=c(outR1, outR2),
                    filter=filterMatch, # Only keep the reads with barcodes
                    compress=TRUE, yieldSize=1e7)
        
}

# Write all the samples to seperate files
mapply(FUN=demultiplexSamples, R1=R1, R2=R2,
       barcode=as.character(barcodes$V1),
       fullMotif=as.character(barcodes$fullMotif))








# clip the barcode of interest from the matched reads
trimLRPatterns(Lpattern=barcode, Rpattern="", )
# Return a fastq files for both pairs with barcode clipped from read











# Setup the simple function for filterFastq
filterFun <- function(x){
        x[match]
}

# Extract reads matching the barcode and output to new fastq file
filterFastq(files=c(fastqPath, pairFile), destinations=outFile,
            filter=filterFun)


