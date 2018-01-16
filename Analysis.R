## Set working directory to source file location!

setwd("~/Dropbox/Studium/Phylogeography Diphasiastrum complanatum")

##——————————————————————————————————————————————————————————————————————————
## Load packages                                             --------------
##——————————————————————————————————————————————————————————————————————————

## Basic stuff
library(stringr) # string manipulation
library(ggplot2)

## Geo
library(maps)
library(mapdata)
library(mapproj)
library(ggmap)

## Genetics
library(Biostrings) # Bio sequence manipulation, here for reading fasta files

## Haplotype Networks
library(stringi) # for generating random strings from set of letters (hack to fool pegas::haplotype)
library(pegas) # Population and Evolutionary Genetics Analysis System. E.g. haploNet()

## Dendrograms
library(phangorn) # upgma clustering, just a wrapper fo hclust(method = "average") (= UPGMA)?
library(ggtree) # plot trees with ggplot2
library(phyloseq) # nice tree plotting with dots
library(hierfstat) # provides Nei's D_A 1983 in function genet.dist()

## Bibliography
library(knitr)
write_bib(x = .packages(), file = "~/Downloads/R-packages.bib", tweak = TRUE, width = NULL, prefix = getOption("knitr.bib.prefix", "R-"))

##——————————————————————————————————————————————————————————————————————————
## Load data                                                 --------------
##——————————————————————————————————————————————————————————————————————————

## save default graphic parameters for later reset with par(.pardefault)
.pardefault <- par(no.readonly = T)

## Won't be necessary everytime:
# source("Alignment cp.R") # execute script preparing the alignment

## loads automatically processed alignment (= Data/Alignment cp processed.fas)
load("Data/Alignment cp processed.RData", verbose = T) # load resulting file, containing
# Metadata # data frame
# cp.alignment.processed

### Replace coordinate column with reconstructed coordinate column
Metadata[,c("Lat", "Lon", "Lat.Lon.Accuracy")] <- Metadata[,c("Lat.constr", "Lon.constr", "Lat.Lon.Accuracy.constr")]

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
setdiff(final.ids, Metadata$Lab.ID)

##————————————————————————————————————————————————————————————————————————————
## Define theme                                          ---------------------
##————————————————————————————————————————————————————————————————————————————

bg.colour <- "#F6F4F0"
line.colour <- "#101010"

pop.colours <- c(Germany = "#F83600", Lithuania = "#006ACD", Kamtchatka = "#0EDD3C", Alaska = "#FFDA00")
pop.9.colours <- rainbow(9)

cols.17.haps <- c("F83600", "006ACD", "7FE700", "FFDA00")

# pop.colors <- c(ByF = "#0000FF", Den = "red", Har = "#003FBF", Lit = "#00BF3F", Sib = "orange", Slo = "#00FF00", ThF = "#007F7F")


##————————————————————————————————————————————————————————————————————————————
## Scope                                                 ---------------------
##————————————————————————————————————————————————————————————————————————————

nrow(MD) # 245 unique specimens (repetitions were excluded in "Alignment cp.R")
table(droplevels(MD$Pop)) # for now 10 populations

nrow(MD[!MD$Is.Outgroup,])
table(droplevels(MD$Pop[!MD$Is.Outgroup])) # for now 10 populations


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
geographic.distances <- distm(coords, coords, fun = distGeo) # distVincentyEllipsoid
# dimnames(geographic.distances) <- list(MD$Pop.landscape, MD$Pop.landscape)
# dimnames(geographic.distances) <- list(MD$Pop, MD$Pop) # can be used for later tree plotting
geographic.dist <- as.dist(geographic.distances, diag = FALSE, upper = FALSE)


## alternatively only the longitudinal distance (along one latitude, averaged)
mean.lats <- aggregate(coords[,"Lat"], by = list(MD$Pop), FUN = mean)
mean.lat <- mean(mean.lats$x) # mean latitude of the mean latitudes of 9 populations:  57.34942
coords.with.fixed.lat <- cbind(coords[,"Lon"], rep(mean.lat, nrow(coords)))
longitudinal.distances <- distm(coords.with.fixed.lat, coords.with.fixed.lat, fun = distGeo) # distVincentyEllipsoid
dimnames(longitudinal.distances) <- list(MD$Pop.landscape, MD$Pop.landscape)
longitudinal.dist <- as.dist(longitudinal.distances, diag = FALSE, upper = FALSE)


##### Cluster pops based on geographical distance
# library(optpart)
# 
# plot(pop.tree <- hclust(d = longitudinal.dist, method = "ward.D2"))
# str(pop.tree)
# tree.pops <- cutree(pop.tree, 5)
# table(tree.pops)


##————————————————————————————————————————————————————————————————————————————
## Evaluate haplotypes from alignment ----------------------------------------
##————————————————————————————————————————————————————————————————————————————

#### Data frame of haplotype descriptors for alignment
## expects: DNAStringSet alignment (better: alignment confined to interesting region)
## returns: Data frame
##
## alternative to stringr: use matchPattern() or vmatchPattern() if you need to find all the occurrences
## (eventually with indels) of a given pattern in a reference sequence or set of sequences.

parse_haplotypes <- function(alignment,
                             snp.regexp = "(?<=GAAAGAAAAA)[AGCT]{1,}(?=AAAA(\\-){1,}CCC)"){
  
  al.strings <- sapply(alignment, FUN = as.character)
  
  snp <- str_extract(al.strings, snp.regexp)
  
  ssr1 <- str_count(str_extract(al.strings, "(?<=TTTC)(AAAT){1,}(?=\\-)"), "AAAT")
  ssr2 <- str_count(str_extract(al.strings, "(?<=\\-)(AT){1,}(?=\\-)"), "AT")
  ssr3 <- str_count(str_extract(al.strings, "(?<=\\-)(ATTT){1,}(?=\\-)"), "ATTT")
  
  haplotype.string <- paste(snp, ssr1, ssr2, ssr3, sep = "_")
  short.haplotype.string <- paste(ssr1, ssr2, ssr3, sep = "-")
  data.frame(SNP = snp, 
             SSR1 = ssr1, 
             SSR2 = ssr2, 
             SSR3 = ssr3, 
             Haplotype = haplotype.string,
             Short.Haplotype = short.haplotype.string,
             File.Names = names(al.strings))
}


#### Compute a step parsimony distance between two haplotype pairs
## expects: two named sample vectors (gt1, gt2) each with names "SNP", "SSR1" … "SSR3"
## returns: distance: count of mutational steps, where A->G counts as 1, and differences in respective SSR repetitions count as mutational steps
ht_distance <- function(ht1, ht2){
  d <- as.numeric(ht1["SNP"] != ht2["SNP"]) # +1 for different SNPs
  ssr.names <- c("SSR1", "SSR2", "SSR3")
  # ht1["SSR1"] <-  ht1["SSR1"] * 2
  # ht2["SSR1"] <- ht2["SSR1"] * 2
  d <- d + sum(abs(ht1[ssr.names] - ht2[ssr.names])) # + sum of absolute SSR count differences
  #d <- d +as.numeric(gt1["Species"] != gt2["Species"])*100
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
## Analysis                                       ----------------------------
##————————————————————————————————————————————————————————————————————————————

#### Do the haplotype reading
## caf alias cp.alignment.final
Haplotypes <- parse_haplotypes(caf)
MD <- cbind(MD, Haplotypes)

#### Write a haplotype table out to csv file
Out.Data <- cbind(Lab.No = MD$Lab.ID,
                  Haplotypes,
                  Pop = MD$Pop,
                  Included.in.big.pops = !MD$Is.Outgroup)
# Out.Data.longer <- Out.Data[match(Metadata.all$Lab.ID, Out.data$Lab.No),]
write.csv(Out.Data, "Data/Haplotypes.csv")
# write.csv(Out.Data.longer, "Data/Haplotypes longer.csv")

#### Adjust the order of levels for plotting
MD$Pop.longitudinal <- factor(MD$Pop.longitudinal, levels = c("Europe.central", "Europe.east", "Siberia.west", "Kamtchatka", "Alaska"))
MD$Pop.landscape <- factor(MD$Pop.landscape, levels = c("Har", "ThF", "ByF", "Slovenia", "Lithuania", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s", "Kamtchatka.central", "Kamtchatka.north", "Alaska"))
MD$Pop.state <- factor(MD$Pop.state, levels = c("Germany", "Slovenia", "Lithuania", "Russia.european", "Ural", "Siberia.west", "Kamtchatka", "Alaska"))
MD$Pop.grouped <- factor(MD$Pop.grouped, levels = c("Europe.central: Har", "Europe.central: ThF", "Europe.central: ByF", "Europe.central: Slovenia", "Europe.east: Lithuania", "Europe.east: St.Petersburg", "Europe.east: Moscow", "Siberia.west: Ural", "Siberia.west: Siberia.west.n", "Siberia.west: Siberia.west.s", "Kamtchatka: Kamtchatka.central", "Kamtchatka: Kamtchatka.north", "Alaska: Alaska"))
MD$Pop <- factor(MD$Pop, levels = c("Germany", "Slovenia", "Lithuania", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s", "Kamtchatka", "Alaska"))

MD$Haplotype <- reorder(MD$Haplotype, MD$SSR1)
MD$Short.Haplotype <- reorder(MD$Short.Haplotype, MD$SSR1)


#### Make haplotype inventory spine plots by region
par(las = 2, mar = c(8, 5, 1, 1))
cls <- function(x) rainbow(length(levels(x)), end = 0.85)

spineplot(Haplotype ~ Pop, data = MD,
          col = cls(MD$Haplotype),
          ylab = "", xlab = "")
labels.pops <- split(MD$Haplotype, MD$Pop) %>%
  sapply(droplevels) %>%
  sapply(levels) %>%
  sapply(rev)
labels.pops

par(.pardefault)

## write textfile with haplotypes for pops
# lapply(labels.pops, write, "pops.txt", append = TRUE, ncolumns = 1000)

table(MD$Haplotype, MD$Pop)

#### build data frame of unique haplotypes
## alternatively, when dealing with multiple species:
## use permutations of Haplotype and Species to return the length of the Lab.ID vector, i. e. the count of the samples
# uhd <- aggregate(Lab.ID ~ Haplotype + Species.mol, data = MD, FUN = length) 

Haplotype.Table <- table(MD$Short.Haplotype, MD$Pop)
Haplotype.Data <- MD[match(rownames(Haplotype.Table), MD$Short.Haplotype), c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Short.Haplotype")]
attr(Haplotype.Table, "class") <- "matrix"
total.count <- rowSums(Haplotype.Table)
Haplotype.Data <- cbind(as.matrix(Haplotype.Table),
                        Haplotype.Data,
                        Total.Count = total.count)
## relevel the Haplotype factor in SSR1 order
Haplotype.Data$Short.Haplotype <- reorder(Haplotype.Data$Short.Haplotype, Haplotype.Data$SSR1)
HD <- Haplotype.Data



## … for narrower ("landscape level") pops: HDL
HT.landscape <- table(MD$Short.Haplotype, MD$Pop.landscape)
HDL <- MD[match(rownames(HT.landscape), MD$Short.Haplotype), c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Short.Haplotype")]
attr(HT.landscape, "class") <- "matrix"
# Total.Count <- rowSums(HT.landscape)
HDL <- cbind(as.matrix(HT.landscape),
             HDL,
             Total.Count = total.count)
## relevel the Haplotype factor in SSR1 order
HDL$Short.Haplotype <- reorder(HDL$Short.Haplotype, HDL$SSR1)

#### specimen distance matrix
## this is quite computation intensive
# specimen.distances <- ht_distance_matrix(MD) # note: there will be NAs in the upper triangle

#### Cluster the specimens
# specimen.dendrogram <- as.dist(specimen.distances) %>%
#   hclust(method = "ward.D2")
# plot(specimen.dendrogram, label = paste(MD$Pop.grouped, MD$Haplotype, sep = " "))

#### Mantel test
# Null hypothesis: these two matrices are uncorrelated
library(ade4) # mantel.rtest(), Data Analysis functions to analyse Ecological and Environmental data in the framework of Euclidean Exploratory methods
mt <- mantel.rtest(as.dist(geographic.dist), as.dist(specimen.distances), nrepet = 9999)
summary(mt)

#### Amova
HD.pops <- HD[, !(names(HD) %in% c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Total.Count", "Short.Haplotype", "Slovenia", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s"))]
amova(samples = HD.pops, distances = as.dist(haplotype.distances), distances = as.dist(haplotype.distances), data.frame(c(1:4)))
 


##————————————————————————————————————————————————————————————————————————————
## Prepare data for Nested clade analysis           --------------------------
##————————————————————————————————————————————————————————————————————————————

# ## Compute a centroid for populations [geosphere::centroid()]
# germany.center <- centroid(as.matrix(MD[MD$Pop == "Germany", c("Lon", "Lat")]))
# lithuania.center <- centroid(as.matrix(MD[MD$Pop == "Lithuania", c("Lon", "Lat")]))
# kamtchatka.center <- centroid(as.matrix(MD[MD$Pop == "Kamtchatka", c("Lon", "Lat")]))
# alaska.center <- centroid(as.matrix(MD[MD$Pop == "Alaska", c("Lon", "Lat")]))
# 
# ## Pop sizes
# table(droplevels(MD$Pop[!MD$Is.Outgroup])) # for now 10 populations
# 
# ## Nexus
# cp.alignment.nca <- cp.alignment.final[!MD$Is.Outgroup]
# names(cp.alignment.nca) <- paste0("DIPcom", MD$Lab.ID, ".", MD$Pop)[!MD$Is.Outgroup]
# writeXStringSet(cp.alignment.nca, filepath = "Data/Alignment NCA.fas")


##————————————————————————————————————————————————————————————————————————————
## Dendrograms plotting      -------------------------------------------------
##————————————————————————————————————————————————————————————————————————————

#### I. Haplotypes
## Subset the haplotype data, e.g. to populations
HD.sub <- HD #[HD$Kamtchatka > 0,]

## Distances
haplotype.distances <- ht_distance_matrix(HD.sub) # note: there will be NAs in the upper triangle
dimnames(haplotype.distances) <- list(HD.sub$Short.Haplotype, HD.sub$Short.Haplotype) # dimnames are passed on as labels to derivated classes

## Distances including outgroup
manual.outgroup <- data.frame(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, SNP = as.factor("A"), SSR1 = 2, SSR2 = 0, SSR3 = 2, Total.count = 2, Haplotype = as.factor("A_2_0_2"), Short.Haplotype = as.factor("2-0-2"))
colnames(manual.outgroup) <- colnames(HD.sub)
HD.out <- rbind(HD.sub, manual.outgroup)
haplotype.distances.out <- ht_distance_matrix(HD.out)
dimnames(haplotype.distances.out) <- list(HD.out$Short.Haplotype, HD.out$Short.Haplotype) # dimnames are passed on as labels to derivated classes

## Make a haplotype abundance matrix based on dropping all other (non abundance and small pop) columns from HD 
pop.matrix <- as.matrix(HD[, !(names(HD) %in% c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Total.Count", "Short.Haplotype", "Slovenia", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s"))]) # "Slovenia", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s"
# rownames(pop.matrix) <- HD$Short.Haplotype # unnecessary, but crucial


## Plot a neighbour joining tree
nj.phylo <- bionj(haplotype.distances.out)
nj.phylo.rooted <- root.phylo(nj.phylo, outgroup = "2-0-2", resolve.root = T)
is.rooted(nj.phylo.rooted)
haplotype.phylo <- phyloseq(otu_table(pop.matrix, taxa_are_rows = T), nj.phylo.rooted)
plot_tree(haplotype.phylo, color = "samples", size = "abundance", label.tips = "taxa_names", sizebase = 2, base.spacing = 0.05) # phyloseq

## Plot a UPGMA tree
upgma.phylo <- upgma(haplotype.distances)
haplotype.phylo <- phyloseq(otu_table(pop.matrix, taxa_are_rows = T), upgma.phylo)
plot_tree(haplotype.phylo, color = "samples", size = "abundance", label.tips = "taxa_names", sizebase = 2, base.spacing = 0.08) +
  scale_colour_manual(values = pop.colours)
 # phyloseq
# cut the tree into 1, 2, 3
# draw x-axis!

## Plot a WARD's tree
ward.phylo <- hclust(as.dist(haplotype.distances), method = "ward.D2")
haplotype.phylo <- phyloseq(otu_table(pop.matrix, taxa_are_rows = T), as.phylo(ward.phylo))
plot_tree(haplotype.phylo, color = "samples", size = "abundance", label.tips = "taxa_names", sizebase = 2, base.spacing = 0.08) +
  scale_colour_manual(values = pop.colours)
# phyloseq

## Note: this is a ggplot object!
# ggtree(haplotype.phylo) +
#   geom_point(aes(x = x + hjust, color = "samples"), na.rm = T) +


#### II. Unrooted tree for populations

## Nei's D_A 1983, (eqn 7) in Takezaki and Nei (1996)
Neis.Data <- MD[!MD$Is.Outgroup, c(Pop = "Pop", # 
                    "SSR1",
                    "SSR2",
                    "SSR3")] # genet.dist expects A data frame containing population of origin as the first column and multi-locus genotypes in following columns
Neis.Data <- droplevels(Neis.Data) # important!
## following step is not necessary, as genet.dist does the factor conversion
Neis.Data[, c("SSR1", "SSR2", "SSR3")] <- lapply(Neis.Data[, c("SSR1", "SSR2", "SSR3")], as.factor)
neis.Da <- genet.dist(Neis.Data, method = "Da", diploid = F)

## Build Neighbour joining dendrogram
## ape::nj
## Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution, 4, 406–425.
## bionj: This function performs the BIONJ algorithm of Gascuel (1997).

neis.Da <- as.matrix(neis.Da)
dimnames(neis.Da) <- list(levels(Neis.Data$Pop), levels(Neis.Data$Pop))
nj.pop.cluster <- bionj(neis.Da) # nj() for classic Saitou amd Nei (1987)
# tiplabels(pop.clust) <- levels(Neis.Data$Pop)
# pop.clust <- hclust(neis.Da, method = "ward.D2")
plot.phylo(nj.pop.cluster, type = "unrooted") # type = "unrooted" # ape::


##————————————————————————————————————————————————————————————————————————————
## Draw map                                              ---------------------
##————————————————————————————————————————————————————————————————————————————



ylims <- c(90, 40) # longitudinal limits for drawing expressed in degrees (range from to)
xlims <- c(NA, NA) # latitudinal limits
# projection.true.degrees <- c(50, 90) # parameters for albers or lambert projection

map.theme <- theme(panel.background = element_rect(fill = "white"),
                   panel.grid.major = element_line(colour = "gray90"))

xscale <- scale_x_continuous(breaks = NULL) # scale ticks for x axis: invisible
Pop <- MD$Pop[!MD$Is.Outgroup]

library(rworldmap)
data(coastsCoarse)
# holarctic <- geom_map(data = coastsCoarse, map = fortify(coastsCoarse), colour = line.colour, fill = "white", map_id = 1:5128) # create a layer of borders, use "worldHires" for publication plotting
holarctic <- borders("world", colour = line.colour, fill = bg.colour) # create a layer of borders, use "worldHires" for publication plotting
specimen.points <- geom_point(aes(x = Lon, y = Lat, group = Pop, color = Pop), size = 2.3, shape = 21, data = droplevels(MD[!MD$Is.Outgroup,])) # MD[!MD$Is.Outgroup,]
map.projection <- coord_map(projection = "stereographic", ylim = ylims)#, parameters = c(60)) # , ylim = ylims)
map <- ggplot() + map.projection + holarctic + specimen.points + map.theme + xscale + scale_colour_manual(values = pop.colours, guide = FALSE)
map

# geom_encircle(aes(x=lon, y=lat), data = places_loc, size = 2, color = "blue")


##————————————————————————————————————————————————————————————————————————————
## Plot haplotype network  ---------------------------------------------------
##————————————————————————————————————————————————————————————————————————————


#### build haplotype network
## as a hacky workaround generate a random alignment with as many different sequences as haplotypes
no.haplotypes <- length(unique(MD$Short.Haplotype))
set.seed(1)
random.strings <- as.list(stri_rand_strings(no.haplotypes, 100, pattern = "[ACTG]"))
random.strings <- as.DNAbin.alignment(ape::as.alignment(random.strings))
names(random.strings) <- HD$Haplotype
ht.d <- ht_distance_matrix(HD)
ht <- pegas::haplotype(random.strings, d = ht.d)
attr(ht, "dimnames")[[1]] <- HD$Short.Haplotype

ht.net <- haploNet(ht, d = ht.d, getProb = T)
attr(ht.net, "freq") <- HD$Total.Count

#### I. Plot the network based on broader (longitudinal) populations
pop.matrix <- as.matrix(HD[, !(names(HD) %in% c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Total.Count", "Short.Haplotype", "Slovenia", "St.Petersburg", "Moscow", "Ural", "Siberia.west.n", "Siberia.west.s"))])

plot(ht.net,
     labels = T,
     threshold = c(0, 0), # no alternative mutation links but smallest distance, 0 otherwise c(1,2)
     size = HD$Total.Count^(1/2)*0.3, # circle sizes
     show.mutation = 1,
     pie = pop.matrix,
     legend = c(10, 10), # coordinates where to draw legend
     bg = pop.colours)

## manual new layout
# layout <- replot()
# replot(layout)


# #### II. Plot the network based on narrower (landsape) populations
# landscape.pop.matrix <- as.matrix(HDL[, !(names(HDL) %in% c("SNP", "SSR1", "SSR2", "SSR3", "Haplotype", "Total.Count", "Short.Haplotype"))])
# c.landscape <- colorRampPalette(c("blue", "green", "yellow", "orange", "red", "purple"))(ncol(landscape.pop.matrix))
# 
# # pop.colors <- c(ByF = "#0000FF", Den = "red", Har = "#003FBF", Lit = "#00BF3F", Sib = "orange", Slo = "#00FF00", ThF = "#007F7F")
# 
# plot(ht.net,
#      labels = T,
#      threshold = c(1,2), # no alternative mutation links but smallest distance, 0 otherwise c(1,2)
#      size = HDL$Total.Count^(1/2)*0.35, # circle sizes
#      show.mutation = 1,
#      pie = landscape.pop.matrix,
#      legend = c(10, 10), # coordinates where to draw legend
#      bg = c.landscape)


##————————————————————————————————————————————————————————————————————————————
## Draw 3D haplotype space                              ---------------------
##————————————————————————————————————————————————————————————————————————————

library(rgl)

plot3d(HD$SSR1, HD$SSR2, HD$SSR3, type = 's', size = total.count^(1/3))


##————————————————————————————————————————————————————————————————————————————
## Model SSR repeats                                      ---------------------
##————————————————————————————————————————————————————————————————————————————
Lon.pos <- MD$Lon
Lon.pos[Lon.pos < 0] <- Lon.pos[Lon.pos < 0] + 360

plot(SSR3 ~ Lon.pos, data = MD)
summary(lm(SSR3 ~ Lon.pos, data = MD))

### SSR1 variation
MD.x.3.4 <- MD[MD$SSR2 == 3 & MD$SSR3 == 4,]
Lon.pos <- MD.x.3.4$Lon
Lon.pos[Lon.pos < 0] <- Lon.pos[Lon.pos < 0] + 360

plot(SSR1 ~ Lon, data = MD.x.3.4)
summary(lm(SSR1 ~ Lon, data = MD.x.3.4))
## nichts geographisches Erkennbar


