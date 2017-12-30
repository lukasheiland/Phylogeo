##——————————————————————————————————————————————————————————————————————————
## Load packages                                              --------------
##——————————————————————————————————————————————————————————————————————————
# library("ape") # for reading FASTA files
library("Biostrings") # guess what
library("sangerseqR") # abif reading, PolyPeakParser,
library("msa") # multiple Sequence alignment


##——————————————————————————————————————————————————————————————————————————
## Functions                                                  --------------
##——————————————————————————————————————————————————————————————————————————

#### Read ab1 files (= abifs)
## expects: full file name of an abif
## returns: named sangerSeq
read_seqs <- function(abif.full){
  seq <- readsangerseq(abif.full) # equivalent to sangerseq(read.abif("file"))
  lab.id <- str_extract(abif.full, "(?<=/)([0-9]{2,3})(?=_)")
  attr(seq, "lab.id") <- lab.id
  return(seq)
}

#### Add gaps to a DNAString to align it to a DNAStringset
## expects: DNAString
## returns: 
align_seq <- function(seq, alignment = cp.alignment.dip.com, reference.index = 1, type = "global-local"){
  pairwiseAlignment(pattern = alignment[reference.index], subject = seq, type = type)
    # type = "global-local" is crucial to align the subject inside the reference pattern
}



##——————————————————————————————————————————————————————————————————————————
## Set "environment" variables                                --------------
##——————————————————————————————————————————————————————————————————————————

source("File management.R") # renames files, and makes sure that the following vectors are usable and valid
cp.files.full <- list.files(cp.seq.directory, full.names = T)
cp.file.ids <- str_extract(cp.files.full, "(?<=/)([0-9]{2,3})(?=_)")


##——————————————————————————————————————————————————————————————————————————
## Load and cleanup original manually edited alignment        --------------
##——————————————————————————————————————————————————————————————————————————

cp.alignment <- readDNAStringSet(file = "Data/Alignment cp original.fas") # package:Biostrings returns DNAStringSet
## alternatively
# cp.alignment <- read.FASTA(file = "Data/Alignment cp manually edited.fas") # package:ape returns DNAbin

to.be.discarded <- grepl("^(cp000|cp999|re|xx)", labels(cp.alignment)) # cp000 are primers etc., cp999 imported from genebank, re repetitions, xx
cp.alignment <- cp.alignment[!to.be.discarded] # drop unwanted sequences


##!! Please make sure that in the alignment prefixes are unique!
is.duplicate <- duplicated(Map("[", strsplit(labels(cp.alignment), "_", fixed = T), 1)) # checks duplicates based on the prefix before "_"
cp.alignment <- cp.alignment[!is.duplicate] ## drop duplicates


##!! attributes are extracted from the labels of form someID_GENspec_labID
##   and added to sequences in the alignment
split.labels <- strsplit(labels(cp.alignment), "_", fixed = T)
attr(cp.alignment, "species") <- sapply(split.labels, function(x) x[2])
attr(cp.alignment, "lab.id") <- sapply(split.labels, function(x) x[3])
## feed length of a split label vector as index, extracting the last of the strings
# attr(cp.alignment, "region") <- sapply(split.labels, function(x) x[length(x)])


## drop everything that is not _DIPcom_ for now
is.dip.com <- attr(cp.alignment, "species") == "DIPcom"
cp.alignment.dip.com <- cp.alignment[is.dip.com] # (vertical subset)
## attributes of the alignment are simple vectors related to the whole aligmment and have to be added anew
split.labels.dip.com <- split.labels[which(is.dip.com)]
attr(cp.alignment.dip.com, "species") <- sapply(split.labels.dip.com, function(x) x[2])
attr(cp.alignment.dip.com, "lab.id") <- dip.com.ids <- sapply(split.labels.dip.com, function(x) x[3])


##——————————————————————————————————————————————————————————————————————————
## Add sequences from abif and scf files to alignment         --------------
##——————————————————————————————————————————————————————————————————————————

## Cascade of logicals to finally obtain two vectors of all new files not yet in the alignment
## fulfilling the criteria:
## 1. edited.
## 2. unedited and there exists none edited with the same lab.id,

are.new.ids <- !(cp.file.ids %in% attr(cp.alignment.dip.com, "lab.id")) # logical vector of whether lab.ids are in the files but not in the manual alignment
are.edited.ids <- str_detect(cp.files, "(?<=_|-)(edit)") # logical vector of whether files are RC'ed and edited 
are.new.edited.ids <- are.new.ids & are.edited.ids
are.new.unedited.ids <- are.new.ids & !are.edited.ids

are.new.unedited.but.ids.exists.in.edited.ids <- are.new.unedited.ids &
  (cp.file.ids %in% intersect(cp.file.ids[are.new.unedited.ids], cp.file.ids[are.new.edited.ids]))
  ## logical vector of whether there exists a file which is (0) new, (2) belongs to the duplicate ids which are in both unedited and edited,
  ## and (1) is TRUE in the logical vector are.new.unedited (to exclude the edited)
are.new.none.edited.ids <-  are.new.unedited.ids & !are.new.unedited.but.ids.exists.in.edited.ids


cp.seqs.none.edited <- cp.files.full[are.new.none.edited.ids] %>% # list of primarySeqs of sangerSeq objects, already RC'ed
  sapply(read_seqs) %>%
  sapply(primarySeq) %>%
  sapply(reverseComplement)

cp.seqs.edited <- cp.files.full[are.new.edited.ids] %>% # list of primarySeqs of sangerSeq objects
  sapply(read_seqs) %>%
  sapply(primarySeq)

cp.seqs <- c(cp.seqs.none.edited, cp.seqs.edited)

cp.seqs.aligned <- sapply(cp.seqs, align_seq) %>% # list of aligned DNAStringsSet
  sapply(subject)  %>% # extracts the subject rather than e.g. the pattern
  sapply(aligned) # extracts the aligned sequence as DNAStringSet
  
  
names(cp.seqs.aligned) <- NULL # necessary for the next command do succeed
cp.seq.alignment <- do.call(c, cp.seqs.aligned)

## construct new names consistent with the manual alignment, so that they are someID_GENspec_labID
added.names <- c(cp.files[are.new.none.edited.ids], cp.files[are.new.edited.ids]) # old file names are expected to consist of labID_etc. anyway
added.ids <- c(cp.file.ids[are.new.none.edited.ids], cp.file.ids[are.new.edited.ids]) # labIDs
new.names <- paste(added.ids, "DIPcom", added.names, sep = "_") # the three strings pasted will yield labID_DIPcom_{labID_includedinoldfilename}

names(cp.seq.alignment) <- new.names

# writePairwiseAlignments((pairwiseAlignment(cp.alignment.dip.com[1], cp.seq.alignment[4])))

cp.alignment.processed <- c(cp.alignment.dip.com, cp.seq.alignment)
# new.alignment <- msa(new.alignment)
# new.alignment <- DNAStringSet(new.alignment)

## visual DNAbin alignment inspection:
# alview(alignment.cp)
# image(alignment.cp)


##——————————————————————————————————————————————————————————————————————————
## Get metadata from master file                               -------------
##——————————————————————————————————————————————————————————————————————————

## Read the csv masterfile
Metadata.raw <- read.csv("Data/Master file.csv",
                         skip = 5,   # skips 5 rows in the .csv file and interprete the 6th as header 
                         check.names = F, # ensures that empty names are not replaced by constructed ones
                         colClasses = c(Lon = "numeric", Lat = "numeric"))
Metadata <- Metadata.raw[, names(Metadata.raw) != ""] # use only rows with set names

##——————————————————————————————————————————————————————————————————————————
## Write data                                                  -------------
##——————————————————————————————————————————————————————————————————————————

writeXStringSet(cp.alignment.processed, "Data/Alignment cp processed.fas")
save(cp.alignment.processed, Metadata, file = "Data/Alignment cp processed.Rdata")
