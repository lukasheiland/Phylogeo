

# haploGen # Simulation of genealogies of haplotypes

##——————————————————————————————————————————————————————————————————————————
## Load packages                                             --------------
##——————————————————————————————————————————————————————————————————————————

## Basic stuff
library(stringr) # string manipulation

## Genetics
library(pegas) # Population and Evolutionary Genetics Analysis System. E.g. haplonet

## Geo

##——————————————————————————————————————————————————————————————————————————
## Load data                                                 --------------
##——————————————————————————————————————————————————————————————————————————

## Won't be necessary everytime:
# source("Alignment cp.R") # execute script preparing the alignment

## loads automatically processed alignment (= Data/Alignment cp processed.fas)
load("Data/Alignment cp processed.RData", verbose = T) # load resulting file, containing
# Metadata # data frame
# cp.alignment.processed

### Replace real coordinates with fake column
Metadata[,c("Lat", "Lon", "Lat.Lon.Accuracy")] <- Metadata[,c("Lat.fake", "Lon.fake", "Lat.Lon.Accuracy.fake")]
### alias Metadata
MD <- Metadata


## in case there are manual edits to the automatically processed alignment,
## get a final alignment file (= Data/Alignment cp final.fas)
cp.alignment.final <- readDNAStringSet(file = "Data/Alignment cp processed final.fas") # package:Biostrings returns DNAStringSet
## make sure that it is sorted by id
final.ids <- sapply(str_split(names(cp.alignment.final), "_"), function(x) x[3])
caf <- cp.alignment.final[order(final.ids)]
ids <- final.ids[order(final.ids)]

## check if ids and Metadata lab.ids are identical
identical(ids, as.character(MD$Lab.ID))


##————————————————————————————————————————————————————————————————————————————
## Scope                                                 ---------------------
##————————————————————————————————————————————————————————————————————————————

nrow(MD) # 252 unique specimens (repetitions were excluded in "Alignment cp.R")
table(droplevels(MD$Pop)) # for now 8 populations


##————————————————————————————————————————————————————————————————————————————
## Evaluate Geodata                                      ---------------------
##————————————————————————————————————————————————————————————————————————————

#### Run some checks first
## check if all specimens have coordinates
coords.present <- !(is.na(MD$Lon) | is.na(MD$Lat))
all(coords.present)

## IDs lacking coordinates:
MD$Lab.ID[!coords.present] # for now: 395 396 397 420 441 455 465 638
MD$Pop.ger[!coords.present] # from everywhere

MD.geo <- MD[coords.present,]
caf.geo <- caf[coords.present]

## check coordinate plausibility
data.frame(MD$Lab.ID, MD$Pop, map.where(database = "world", MD$Lon, MD$Lat))
## Note that 203 in "ByF" is actually in Czech republic

#### Geodistance
library(geosphere)

## a great circle distance was computed between each sample by means of the Vincenty ellipsoid distance method provided by the {geosphere} package.
## default unit is meters
coords <- cbind(Lon = MD$Lon, Lat = MD$Lat)
geographic.distances <- distm(coords, coords, fun = distVincentyEllipsoid)
dimnames(geographic.distances) <- list(MD$Pop.landscape, MD$Pop.landscape)
# dimnames(geographic.distances) <- list(MD$Pop, MD$Pop) # can be used for later tree plotting
geographic.dist <- as.dist(geographic.distances, diag = FALSE, upper = FALSE)


## alternatively only the longitudinal distance (along one latitude, averaged)
mean.lats <- aggregate(coords[,"Lat"], by = list(MD$Pop.longitudinal), FUN = mean)
mean.lat <- mean(mean.lats$x) # mean latitude of the mean latitudes of 9 populations:  57.34942
coords.with.fixed.lat <- cbind(coords[,"Lon"], rep(mean.lat, nrow(coords)))
longitudinal.distances <- distm(coords, coords, fun = distVincentyEllipsoid)
dimnames(longitudinal.distances) <- list(MD$Pop.landscape, MD$Pop.landscape)
longitudinal.dist <- as.dist(longitudinal.distances, diag = FALSE, upper = FALSE)

##### Cluster pops
library(optpart)


plot(pop.tree <- hclust(d = longitudinal.dist, method = "ward.D2"))
str(pop.tree)
tree.pops <- cutree(pop.tree, 5)
table(tree.pops)


##————————————————————————————————————————————————————————————————————————————
## Evaluate haplotypes from alignment ----------------------------------------
##————————————————————————————————————————————————————————————————————————————

# Use matchPattern or vmatchPattern if you need to find all the occurrences
#  (eventually with indels) of a given pattern in a reference sequence or set of sequences.


#### Data frame of haplotype descriptors for alignment
## expects: DNAStringSet alignment (better: alignment confined to interesting region)
## returns: Data frame

parse_haplotypes <- function(alignment,
                             snp.regexp = "(?<=GAAAGAAAAA)[AGCT]{1,}(?=AAAA(\\-){1,}CCC)"){
  
  al.strings <- sapply(alignment, FUN = as.character)
  
  snp <- str_extract(al.strings, snp.regexp)
  
  ssr1 <- str_count(str_extract(al.strings, "(?<=TTTC)(AAAT){1,}(?=\\-)"), "AAAT")
  ssr2 <- str_count(str_extract(al.strings, "(?<=\\-)(AT){1,}(?=\\-)"), "AT")
  ssr3 <- str_count(str_extract(al.strings, "(?<=\\-)(ATTT){1,}(?=\\-)"), "ATTT")
  
  haplotype.string <- paste(snp, ssr1, ssr2, ssr3, sep = "_")
  data.frame(SNP = snp, 
             SSR1 = ssr1, 
             SSR2 = ssr2, 
             SSR3 = ssr3, 
             Haplotype = haplotype.string,
             File.Names = names(al.strings))
}


#### Compute a simple step distance between two haplotype pairs
## expects: two named sample vectors (gt1, gt2) each with names "SNP", "SSR1" … "SSR3"
## returns: distance: count of mutational steps, where A->G counts as 1, and differences in respective SSR repetitions count as mutational steps
ht_distance <- function(ht1, ht2){
  d <- as.numeric(ht1["SNP"] != ht2["SNP"]) # +1 for different SNPs
  ssr.names <- c("SSR1", "SSR2", "SSR3")
  d <- d + sum(abs(ht1[ssr.names] - ht2[ssr.names])) # + sum of absolute SSR count differences
  #d <- d +as.numeric(gt1["Species"] != gt2["Species"])*100
  d
}

#### Compute a euclidean distance between two haplotype pairs
## expects: two named sample vectors (gt1, gt2) each with names "SNP", "SSR1" … "SSR3"
## returns: distance: count of mutational steps, where A->G counts as 1, and differences in respective SSR repetitions count as mutational steps
euclid_ht_distance <- function(ht1, ht2){
  snp.d <- as.numeric(ht1["SNP"] != ht2["SNP"]) # +1 for different SNPs
  ssr1.d <- ht1["SSR1"] -  ht2["SSR1"]
  ssr2.d <- ht1["SSR2"] -  ht2["SSR2"]
  ssr3.d <- ht1["SSR3"] -  ht2["SSR3"]
  
  d <- sqrt(snp.d^2 + ssr1.d^2 + ssr2.d^2 + ssr3.d^2) # + sum of squared SSR count differences
  d
}

#### Build a distance matrix for a dataset
## expects: a data.frame with columns: "UID", "SNP", "SSR1" … "SSR3"
## returns: a lower triangle distance matrix
ht_distance_matrix <- function(df, fun = "ht_distance"){
  ## No. of types in data.frame
  n <- nrow(df)
  ## this simply produces a 2 row matrix of possible unique combinations for n values, where n = nrow(dataframe)
  ## combn() achieves "triangular" comparison: of x values, n elements chosen
  comb.pairs <- combn(x = nrow(df), 2)
  ## compute distances
  distances <- apply(comb.pairs,
                     MARGIN = 2,
                     FUN = function(x)  ht_distance(df[x[1],], df[x[2],]))
  ## create empty matrix of n by n, set dimnames to the unique ID
  dist.matrix <- matrix(rep(NA, n^2), nrow = n, dimnames = list(df$Lab.ID, df$Lab.ID))
  ## fill lower triangle of distance matrix
  dist.matrix[lower.tri(dist.matrix)] <- distances
  dist.matrix
}

#### Build a distance matrix for a dataset
## expects: a data.frame with columns: "UID", "SNP", "SSR1" … "SSR3"
## returns: a lower triangle distance matrix
euclid_ht_distance_matrix <- function(df){
  ## No. of types in data.frame
  n <- nrow(df)
  ## this simply produces a 2 row matrix of possible unique combinations for n values, where n = nrow(dataframe)
  ## combn() achieves "triangular" comparison: of x values, n elements chosen
  comb.pairs <- combn(x = nrow(df), 2)
  ## compute distances
  distances <- apply(comb.pairs,
                     MARGIN = 2,
                     FUN = function(x)  euclid_ht_distance(df[x[1],], df[x[2],]))
  ## create empty matrix of n by n, set dimnames to the unique ID
  dist.matrix <- matrix(rep(NA, n^2), nrow = n, dimnames = list(df$Lab.ID, df$Lab.ID))
  ## fill lower triangle of distance matrix
  dist.matrix[lower.tri(dist.matrix)] <- distances
  dist.matrix
}


##————————————————————————————————————————————————————————————————————————————
## Analysis                                       ----------------------------
##————————————————————————————————————————————————————————————————————————————

#### do the haplotype reading
## caf alias cp.alignment.final
Haplotypes <- parse_haplotypes(caf)
MD <- cbind(MD, Haplotypes)

#### adjust the order of levels for plotting
MD$Pop.longitudinal <- factor(MD$Pop.longitudinal,
                                    levels = levels(MD$Pop.longitudinal)[c(1, 2, 4, 5, 6, 3)])
MD$Haplotype <- reorder(MD$Haplotype,
                              MD$SSR1)
MD$Pop.landscape <- factor(MD$Pop.landscape,
                                 levels = levels(MD$Pop.landscape)[c(1, 3, 9, 2, 8, 5, 6, 10, 7, 4)])
MD$Pop.state <- factor(MD$Pop.state,
                                 levels = levels(MD$Pop.state)[c(1, 2, 7, 4, 5, 8, 6, 3)])

#### Make haplotype inventory spine plots by region
par(las = 2, mar = c(8, 5, 1, 1))
cls <- function(x) rainbow(length(levels(x)), end = 0.85)


spineplot(droplevels(Haplotype) ~ Pop.landscape, data = MD,
          col = cls(MD$Haplotype),
          ylab = "", xlab = "")
labels.landscape <- split(MD$Haplotype, MD$Pop.landscape) %>%
  sapply(droplevels) %>%
  sapply(levels) %>%
  sapply(rev)
lapply(labels.landscape, write, "landscape.txt", append = TRUE, ncolumns = 1000)


spineplot(Haplotype ~ Pop.state, data = MD,
          col = cls(MD$Haplotype),
          ylab = "", xlab = "")
labels.state <- split(MD$Haplotype, MD$Pop.state) %>%
  sapply(droplevels) %>%
  sapply(levels) %>%
  sapply(rev)
lapply(labels.state, write, "state.txt", append = TRUE, ncolumns = 1000)


spineplot(Haplotype ~ Pop.longitudinal, data = MD,
          col = cls(MD$Haplotype),
          ylab = "", xlab = "")
labels.longitudinal <- split(MD$Haplotype, MD$Pop.longitudinal) %>%
  sapply(droplevels) %>%
  sapply(levels) %>%
  sapply(rev)
lapply(labels.longitudinal, write, "longitudinal.txt", append = TRUE, ncolumns = 1000)

table(MD$Haplotype, MD$Pop.longitudinal)

#### Cluster the distance!
ht.distances <- ht_distance_matrix(MD) # note: there will be NAs in the upper triangle
ht.dendrogram <- as.dist(ht.distances) %>%
  hclust(method = "ward.D2")
plot(ht.dendrogram, label = paste(MD$Pop.landscape, MD$Haplotype, sep = " "))


euclid.ht.distances <- euclid_ht_distance_matrix(MD)
euclid.ht.dendrogram <- as.dist(euclid.ht.distances) %>%
  hclust(method = "ward.D2")
plot(euclid.ht.dendrogram, label = paste(MD$Pop.landscape, MD$Haplotype, sep = " "))



##————————————————————————————————————————————————————————————————————————————
## Map drawing                                           ---------------------
##————————————————————————————————————————————————————————————————————————————
library(maps)
library(mapdata)
library(mapproj)
library(ggmap)


ylims <- c(40, 90) # longitudinal limits for drawing expressed in degrees (range from to)
xlims <- c(NA, NA) # latitudinal limits
projection.true.degrees <- c(50, 90)


map.theme <- theme(panel.background = element_rect(fill = "black"),
                   panel.grid.major = element_line(colour = "gray20"))

xscale <- scale_x_continuous(breaks = NULL)
Pop <- MD$Pop.longitudinal

holarctic <- borders("world", colour = "gray85", fill = "gray93") # create a layer of borders, use "worldHires" for publication plotting
# rivers <- borders("rivers", colour = "gray80", fill = "black") # create a layer of borders, use "worldHires" for publication plotting
specimen.points <- geom_point(aes(x = Lon, y = Lat, group = Pop, color = Pop), size = 2, data = MD)
map.projection <- coord_map(projection = "stereographic", ylim = ylims)#, parameters = c(60)) # , ylim = ylims)
map <- ggplot() + map.projection + holarctic + specimen.points + map.theme + xscale
map



# ## Attributes to plot parameters
# vertex.sizes <- D$Count*3 # set as Area
# edge.widths <- replace(E(D.graph)$distance, E(D.graph)$distance == 0, 0.3)
# # edge.widths <- 1.3/edge.widths
# edge.widths <- 1/edge.widths

## plot haplo network on map
#plot(D.graph, layout=jitter(coord, 5), add = TRUE, rescale = F,
#     vertex.shape = "circle", # "none" und vertex.label.dist = 0 für nur Test
#     vertex.size=vertex.sizes,
#     vertex.frame.color = "black",
#     vertex.color = adjustcolor("green2", alpha.f = 0.3),
#     vertex.label.color="black",
#     vertex.label.cex=0.5,
#     vertex.label.dist=vertex.sizes*0.02+0.2,
#     vertex.label.degree= -pi/2, # default
#     edge.arrow.mode=0,
#     edge.label.cex=0.6,
#     edge.label.color="black",
#     edge.color= adjustcolor("green", alpha.f = 0.3),
#     edge.width=edge.widths,
#     edge.label=E(D.graph)$distance,
#     edge.curved=0.2,
#     edge.loop.angle=pi/2)



########################################

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











##————————————————————————————————————————————————————————————————————————————
## OLD: Define distance functions                             ----------------
##————————————————————————————————————————————————————————————————————————————

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

##————————————————————————————————————————————————————————————————————————————
## Data inspection -----------------------------------------------------------
##————————————————————————————————————————————————————————————————————————————


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
zehn57 <- subset(MD, Haplotype == "A_10_5_7")
table(zehn57$Pop.landscape)


##————————————————————————————————————————————————————————————————————————————
## Dendrogram plotting -------------------------------------------------------
##————————————————————————————————————————————————————————————————————————————
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




##————————————————————————————————————————————————————————————————————————————
## Confine alignment to relevant region  -------------------------------------
##————————————————————————————————————————————————————————————————————————————

# probably for haplonet


##————————————————————————————————————————————————————————————————————————————
## Plot haplotype networks ---------------------------------------------------
##————————————————————————————————————————————————————————————————————————————


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












##————————————————————————————————————————————————————————————————————————————
## Get metadata from haplotypes from align                      --------------
##————————————————————————————————————————————————————————————————————————————

#### "meta data frame" from a haplotype object
#### extracting Region etc. from associated alignment
get.meta.df <- function(haplotypes, snp.pos = 1){
  al <- get(attr(haplotypes, "from")) # get alignment the haplotypes were derived from
  
  ## extract first indices of each haplotype in alignment
  hi <- sapply(attr(haplotypes, "index"), "[", 1)
  al.h <- al[hi, ] # alignment comprising only unique haplotypes
  species.h <- attr(al, "species")[hi] # this is the way to get alignment attributes in the metadata
  MD <- cbind(hap.factors(al.h), Species = species.h)
  rownames(MD) <- attr(haplotypes, "dimnames")[[1]] # this is problematic in case there are different haplotypes with identical haplotype descriptor, e.g. other SNPs
  MD
}
