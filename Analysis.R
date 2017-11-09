##-----------------------------------------------------------------------------
## Data preparation
##-----------------------------------------------------------------------------

source("Preparation.R") # execute file with data preparation

load("Data/Alignment_processed.RData", verbose = T) # load resulting file, containing some data frames
## EData.RData for data extracted from Excel file
## where DC is an alias for URHD.dc

##-----------------------------------------------------------------------------
## Define distance functions
##-----------------------------------------------------------------------------

#### Compute a distance between two haplotype pairs
## expects: two named sample vectors with names "SNP", "SSR1" … "SSR3"
## returns: distance: count of mutational steps, where A->G counts as 1, and differences in respective SSR repetitions count as mutational steps
distance <- function(gt1, gt2){
  d <- as.numeric(gt1["SNP"] != gt2["SNP"]) # +1 for different SNPs
  ssr.names <- c("SSR1", "SSR2", "SSR3")
  d <- d + sum(abs(gt1[ssr.names] - gt2[ssr.names])) # + sum of absolute SSR count differences
  #d <- d +as.numeric(gt1["Species"] != gt2["Species"])*100
  d
}

#### Build a distance matrix for a dataset
## expects: a data.frame with columns: "UID", "SNP", "SSR1" … "SSR3"
## returns: a lower triangle distance matrix
distance.matrix <- function(df){
  ## No. of types in data.frame
  n <- nrow(df)
  ## this simply produces a 2 row matrix of possible unique combinations for n values, where n = nrow(UHD)
  ## combn() achieves "triangular" comparison: of x values, n elements chosen
  comb.pairs <- combn(x = nrow(df), 2)
  ## compute distances
  distances <- apply(comb.pairs, MARGIN = 2, FUN = function(x) distance(df[x[1],], df[x[2],]))
  ## create empty matrix of n by n, set dimnames to the unique ID
  dist.matrix <- matrix(rep(NA, n^2), nrow = n, dimnames = list(df$UID, df$UID))
  ## fill lower triangle of distance matrix
  dist.matrix[lower.tri(dist.matrix)] <- distances
  dist.matrix
}

#### Distances can also be strucured in a data.frame
## expects: a data.frame with columns: "UID", "SNP", "SSR1" … "SSR3"
## returns: a data.frame of paired haplotypes in the first two rows and respective distances
distance.df <- function(df){
  ## No. of types in data.frame
  n <- nrow(df)
  ## this simply produces a 2 row matrix of possible unique combinations for n values, where n = nrow(UHD)
  ## combn() achieves "triangular" comparison: of x values, n elements chosen
  comb.pairs <- combn(x = nrow(df), 2)
  distances <- apply(comb.pairs, MARGIN = 2, FUN = function(x) distance(df[x[1],], df[x[2],]))
  DDF <- data.frame(df$UID[comb.pairs[1,]], df$UID[comb.pairs[2,]], distances)
  names(DDF) <- c("From", "To", "Distances")
  DDF
}

##### invoke some distance matrices
uhd.distances <- distance.matrix(UHD)
urhd.distances <- distance.matrix(URHD)
dc.distances <- distance.matrix(DC)
dc.ger.distances <- distance.matrix(DC.ger)

##-----------------------------------------------------------------------------
## Data inspection
##-----------------------------------------------------------------------------
germany <- c("Har", "ByF", "ThF", "Alp", "LHe", "OdW")

table(HD$Region)

## Total
nrow(HD) # total no. of records
length(levels(HD$Species)) # total no. of taxa
length(levels(HD$Genotype)) # total no. of genotypes
droplevels(HD$Genotype[HD$Region %in% germany]) # haplotypes in Germany: 14
table(HD$Region) # by region
table(HD$Region) # by region

## DC (DIPcom) Data
DC[c("Region", "Count")]
sum(DC$Count) # no. of records
sum(DC$Count[DC$Region %in% germany])
levels(DC$Genotype) # haplotypes
levels(DC$Genotype[DC$Region %in% germany]) # haplotypes in Germany


## Biogeography
table(DC$Region) # Uniques by region

recordsbyregion <- c(aggregate(DC$Count, by = list(DC$Region), sum)[[2]])
table(DC$Region)/recordsbyregion  # Uniques per samples by region

## count of unique haplotypes in a region are probably proportional to record count
plot(as.numeric(table(DC$Region)) ~ recordsbyregion)
summary(lm(as.numeric(table(DC$Region)) ~ recordsbyregion)) # **, R squared = 0.8

## Average distance between regions
# ddf <- distance.df(DC)
# ddf$From <- DC$Region[match(ddf$From, DC$UID)]
# ddf$To <- DC$Region[match(ddf$To, DC$UID)]
# regional.distances <- aggregate(Distances ~ From + To, data = ddf, mean)

#### Main Haplotypes

vierzehn34 <- subset(HD, Genotype == "A_14_3_4")
zehn57 <- subset(HD, Genotype == "A_10_5_7")

table(zehn57$Region, zehn57$Species)
table(vierzehn34$Region, vierzehn34$Species)


##-----------------------------------------------------------------------------
## Dendrogram plotting
##-----------------------------------------------------------------------------
library("dendextend") # for customizing dendrograms

## plot dendrograms of hierarchichal, agglomerative cluster analyses from distance matrix
plot(hclust(as.dist(uhd.distances)))
plot(hclust(as.dist(urhd.distances)), label = URHD$Region)

#### DIPcom ~ Region
dend <- hclust(as.dist(dc.distances), method = "ward.D2") # best truthful method for building nice groups!
dend <- as.dendrogram(dend)

p.sizes <- (0.8*DC$Count^(1/2))[order.dendrogram(dend)]
p.colors <- DC$Region[order.dendrogram(dend)]

rc <- c(ByF = "#0000FF", Den = "red", Har = "#003FBF", Lit = "#00BF3F", Sib = "orange", Slo = "#00FF00", ThF = "#007F7F")

par(mar=c(2, 1, 2, 9))


dend %>%
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", p.sizes) %>%  # node point size
  # set("labels_col", DC$Region) %>% # change color
  # set("labels_col", c("black", "red", "grey"), k=3) %>% # color labels by specifying the number of cluster (k)
  set("leaves_col", rc[p.colors]) %>% # node point color
  hang.dendrogram() %>%
  plot(main = "Diphasiastrum complanatum", horiz = T)

rect.dendrogram(dend, k=3, border = 8, lty = 5, lwd = 1, horiz = T)


##-----------------------------------------------------------------------------
## Evaluate haplotypes from alignment
##-----------------------------------------------------------------------------

library(pegas)
library(stringr)

#### Subset an alignment
## expects: an alignment (DNAbin) and subsetting parameters
## returns: a (vertically and/or horizontally) subset alignment
subset.alignment <- function(alignment, from = 1, to = ncol(alignment), species = F){
  is.v.subset <- attr(alignment, "species") == species | (species == F)
  al <- alignment[is.v.subset, ]
  al <- as.matrix(al)[, from:to]
  attr(al, "lab.no") <- attr(alignment, "lab.no")[is.v.subset]
  attr(al, "species") <- attr(alignment, "species")[is.v.subset]
  attr(al, "region") <- attr(alignment, "region")[is.v.subset]
  al
}

#### df of haplotype descriptors for alignment
## expects: alignment (better: alignment confined to interesting region)
## returns: df
hap.factors <- function(al, snp.pos = 1){
  snp <- as.character(al[,snp.pos])
  snp <- toupper(snp)
  
  al.strings <- sapply(al, FUN = function(x) paste(as.character(x), collapse=""))
  ssr1 <- str_count(str_extract(al.strings, "(?<=tttc)(aaat){1,}(?=\\-)"), "aaat")
  ssr2 <- str_count(str_extract(al.strings, "(?<=\\-)(at){1,}(?=\\-)"), "at")
  ssr3 <- str_count(str_extract(al.strings, "(?<=\\-)(attt){1,}(?=\\-)"), "attt")
  
  haplotype.string <- paste(snp, ssr1, ssr2, ssr3, sep = "_")
  data.frame(SNP = snp,
             SSR1 = ssr1,
             SSR2 = ssr2,
             SSR3 = ssr3,
             Haplotype = haplotype.string)
}

#### "meta data frame" from a haplotype object
#### extracting Region etc. from associated alignment
meta.df <- function(haplotypes, snp.pos = 1){
  al <- get(attr(haplotypes, "from")) # get alignment the haplotypes were derived from

  ## extract first indices of each haplotype in alignment
  hi <- sapply(attr(haplotypes, "index"), "[", 1)
  al.h <- al[hi, ] # alignment comprising only unique haplotypes
  species.h <- attr(al, "species")[hi] # this is the way to get alignment attributes in the metadata
  MD <- cbind(hap.factors(al.h), Species = species.h)
  rownames(MD) <- attr(haplotypes, "dimnames")[[1]] # this is problematic in case there are different haplotypes with identical haplotype descriptor, e.g. other SNPs
  MD
}

####
##
##
region.matrix <- function(haplotypes){
  al <- get(attr(haplotypes, "from")) # get alignment the haplotypes were derived from
  hr.table <- table(hap.factors(al)$Haplotype, attr(al, "region"))
  hr.table[meta.df(haplotypes)$Haplotype,]
}

####
##
##
species.matrix <- function(haplotypes){
  al <- get(attr(haplotypes, "from")) # get alignment the haplotypes were derived from
  hs.table <- table(hap.factors(al)$Haplotype, attr(al, "species"))
  hs.table[meta.df(haplotypes)$Haplotype,]
}

#### Returns a distance matrix for a haplotype {pegas} object
## expects: haplotype object
## returns: lower triangle distance matrix
haplotype.distances <- function(haplotypes){

}


##-----------------------------------------------------------------------------
## Plot haplotype networks
##-----------------------------------------------------------------------------

#### Layout from NMDS!

#### Plot DIPcom ~ Region
al <- subset.alignment(alignment, from = 364, to = 543, species = "DIPcom") # subset to homogeneous and relevant part of the alignment
h <- haplotype(al)
MD <- meta.df(h) ## why again?
attr(h, "dimnames")[[1]] <- MD$Haplotype
rm <- region.matrix(h)
dm <- distance.matrix(MD)

## Plotting parameters
# ec <- colorRampPalette(c("blue", "green"))(5)
region.colors <- c(ByF = "#0000FF", Den = "red", Har = "#003FBF", Lit = "#00BF3F", Sib = "orange", Slo = "#00FF00", ThF = "#007F7F")

hnet <- haploNet(h, d = dm)
plot(hnet,
     threshold = 0, # no alternative mutation links
     size = attr(hnet, "freq")^(1/2)*0.6,
     pie = rm,
     legend = c(-10, 5), # coordinates
     bg = region.colors)
layout <- replot()
replot(layout)
####--------------
#### Plot DIPalp ~ Region
al <- subset.alignment(alignment, from = 364, to = 543, species = "DIPalp") # subset to homogeneous and relevant part of the alignment
h <- haplotype(al)
MD <- meta.df(h)
attr(h, "dimnames")[[1]] <- MD$Haplotype
rm <- region.matrix(h)
dm <- distance.matrix(MD)

## Plotting parameters
# ec <- colorRampPalette(c("blue", "green"))(5)
region.colors <- c(ByF = "#0000FF", Den = "red", Har = "#003FBF", Lit = "#00BF3F", Sib = "orange", Slo = "#00FF00", ThF = "#007F7F")

hnet <- haploNet(h, d = dm)
plot(hnet,
     threshold = 0, # no alternative mutation links
     size = attr(hnet, "freq")^(1/2)*0.6,
     pie = rm,
     legend = T, # coordinates
     bg = region.colors)


#### Plot All Haplotypes ~ Species
al.all <- subset.alignment(alignment, from = 364, to = 543, species = F)
h.all <- haplotype(al.all)
MD.all <- meta.df(h.all)
attr(h.all, "dimnames")[[1]] <- MD.all$Haplotype
dm.all <- distance.matrix(MD.all)
sm <- species.matrix(h.all)

sp.colors <- c(DIPalp = "yellow", DIPcom = "skyblue", DIPiss = "green", DIPoel = "orange", DIPtri = "red", DIPzei = "violet")

hnet.all <- haploNet(h.all, d = dm.all)
plot(hnet.all,
     #threshold = 0, # no alternative mutation links
     size = attr(hnet.all, "freq")^(1/2),
     legend = T,
     pie = sm,
     bg = sp.colors,
     scale.ratio = 1.6)