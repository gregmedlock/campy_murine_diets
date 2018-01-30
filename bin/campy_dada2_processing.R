library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)

path <- "~/Documents/projects/campy_murine_diets/data/20170224_16S-35570571/"
fns <- list.files(path)
fns

### Load forward and reverse reads
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# NOTE: reads are 250bp instead of 300 this time
plotQualityProfile(fnFs[[3]])
# Looks like first 10 bases and last 10 should be trimmed
plotQualityProfile(fnRs[[3]])
# reverse reads should be trimmed before first 10 and after 150

filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(240,160), trimLeft=c(10,10), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}
# How does looser maxEE affect output? (2,2) was old


### Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

derepFs[[1]]

# learn error rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/projects/campy_murine_diets/data/silva_nr_v123_train_set.fa.gz")
unname(head(taxa))

# create and save phyloseq object
library(phyloseq); packageVersion("phyloseq")

samples.out <- rownames(seqtab.nochim)
plate <- sapply(strsplit(samples.out, "-"), `[`, 1)
well <- sapply(strsplit(samples.out, "-"), `[`, 2)
samdf <- data.frame(Plate=plate, Well=well)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

# Save ps object
saveRDS(ps,"~/Documents/projects/campy_murine_diets/data/campy_phyloseq_obj_less_strict.rds")
# To read later, do:
#ps <- readRDS("path/to/ps.rds")
