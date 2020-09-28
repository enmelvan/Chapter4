# ensures repeatability
set.seed(123) 

# phylum-level analysis
library(beemStatic)
library(ggraph)
library(phyloseq)
library(microbiome)
library(knitr)


# Drew: coding based on relative file paths saves time
my.path       <- '/Users/drew/Desktop/ena_chapter3/Chapter3-master/'
path_for_otu  <- paste0(my.path,'data/Turnbaugh/otu_table.csv')
path_for_taxa <- paste0(my.path,'data/Turnbaugh/taxonomy.csv')
path_for_meta <- paste0(my.path,'data/Turnbaugh/metadata.csv')


#GROUP BY TAXONOMY LEVELS
#LOAD YOUR DATA
otu.file <- read.csv(path_for_otu, row.names = 1)
tax.file <- read.csv(path_for_taxa, row.names = 1)
meta.file <- read.csv(path_for_meta, row.names = 1)


#COVERT IT TO PHYLOSEQ
my_data = phyloseq(
     otu_table(as.matrix(otu.file), taxa_are_rows=T),
     tax_table(as.matrix(tax.file)),
     sample_data(meta.file)
   )


#GROUP BY TAXONOMY LEVEL (INPUT HAS TO BE PHYLOSEQ)
my_phylum <- tax_glom(my_data, taxrank="Phylum")
my_class <- tax_glom(my_data, taxrank="Class")
my_order <- tax_glom(my_data, taxrank="Order")
my_family <- tax_glom(my_data, taxrank = "Family")
my_genus <- tax_glom(my_data, taxrank = "Genus")
my_species <- tax_glom(my_data, taxrank = "Species")


# Drew: I think the key issue is the number of samples with the taxon present, 
# not the mean abundance, so I calculate number of samples with each OTU.
# There are four frequently observed phyla: Otu0001, Otu0044, Otu0076, Otu0339.
otu.nsample <- function(phyloseq_object) {
  n <- c()
  for(i in 1:nrow(phyloseq_object@otu_table@.Data)) {
    temp <- phyloseq_object@otu_table@.Data[i,]
    n[i] <- length(temp[temp>0])
  }
  names(n) <- row.names(phyloseq_object@otu_table@.Data)
  return(n)
}
otu.nsample(my_phylum)
names(otu.nsample(my_phylum))[otu.nsample(my_phylum)>=142/2]


# Drew: I pruned all but the 4 most frequently observed taxa
taxon_sub <- prune_taxa(names(otu.nsample(my_phylum))[otu.nsample(my_phylum)>=142/2],my_phylum)


# Drew: You should add poisson noise to counts, not proportions. I created a new
# object for the noisy data (taxon_sub_rand) and updated the counts directly, 
# with lambda equal to the observed number of counts in the sample.
taxon_sub_rand <- taxon_sub
for(i in 1:ncol(taxon_sub_rand@otu_table@.Data))
  taxon_sub_rand@otu_table@.Data[,i] <- rpois(rep(1,nrow(taxon_sub_rand@otu_table@.Data)),
                                         taxon_sub_rand@otu_table@.Data[,i])


#TRANSFORM TO RELATIVE ABUNDANCE. Drew: I didn't change your code here, but 
# since I am taking this step after pruning taxa and creating random data, the 
# function make_relative() in the library funrar is no longer needed. Everything
# sums to 1 already.
transform_to_rel_ab <- function(phyloseq_object) {
  rel_ab <- transform_sample_counts(phyloseq_object, function(x) x/sum(x))
  return(rel_ab)
}
taxon_rel      <- transform_to_rel_ab(taxon_sub)
taxon_rel_rand <- transform_to_rel_ab(taxon_sub_rand)


# Drew: I never needed to uese these convenience functions. I just extracted the
# abundance data directly below.
# EXTRACT ABUNDANCE MATRIX FROM PHYLOSEQ OBJECT
#extract_abundance <- function(phyloseq_object){
#  OTU2 = as(otu_table(phyloseq_object), "matrix")
#  # transpose if necessary
#  if(taxa_are_rows(phyloseq_object)){OTU2 <- t(OTU2)}
#  # Coerce to data.frame
#  OTUdf2 = as.data.frame(OTU2)
#  return(OTUdf2)
#}
#taxonr_data  <- extract_abundance(taxon_rel)
#taxonr_rand  <- extract_abundance(taxon_rel_rand)


#Lean analysis (seems to work fine, but full convergence not achieved after 1000 iterations)
lean_names    <- row.names(meta.file)[meta.file$Var == "Lean"]
mySubset      <- taxon_rel@otu_table@.Data[,colnames(taxon_rel@otu_table@.Data) %in% lean_names]
mySubset_rand <- taxon_rel_rand@otu_table@.Data[,colnames(taxon_rel_rand@otu_table@.Data) %in% lean_names]
lean.fit <- func.EM(mySubset,mySubset_rand,max.iter = 1e3,epsilon=1e-2)
diagnoseFit(lean.fit, mySubset, annotate = TRUE) # R^2 > 0.996
showInteraction(lean.fit,mySubset)
diagnoseBiomass(lean.fit)                             # convergence looks good
hist(beem2biomass(lean.fit))

#obese analysis (seems to work fine, but full convergence not achieved after 1000 iterations)
obese_names     <- row.names(meta.file)[meta.file$Var == "Obese"]
mySubset      <- taxon_rel@otu_table@.Data[,colnames(taxon_rel@otu_table@.Data) %in% obese_names]
mySubset_rand <- taxon_rel_rand@otu_table@.Data[,colnames(taxon_rel_rand@otu_table@.Data) %in% obese_names]
obese.fit <- func.EM(mySubset,mySubset_rand,max.iter = 1e3,epsilon=1e-2)
diagnoseFit(obese.fit, mySubset, annotate = TRUE) # R^2 > 0.996
showInteraction(obese.fit,mySubset)
diagnoseBiomass2(obese.fit)                             # convergence looks good
hist(beem2biomass(obese.fit))


# modified diagnose biomass function
diagnoseBiomass2 <- function (beem.out, true.biomass = NA, alpha = 0.1, ...) 
{
  trace.m <- beem.out$trace.m
  nIter <- ncol(trace.m)
  if (any(is.na(true.biomass))) {
    true.m <- trace.m[, nIter]
  }
  else {
    true.m <- true.biomass
  }
  rel.err.m <- (t((trace.m - true.m)/true.m)) * 100
  col <- rep(rgb(0, 0, 0, alpha), nrow(trace.m))
  col[beem.out$sample2rm] <- rgb(1, 1, 0, alpha)
  matplot(rel.err.m[floor(nIter/2):nIter], type = "l", xlab = "Iterations", ylab = "Relative difference (%)", 
          main = "Biomass trace", lty = 1, lwd = 3, col = col, 
          ...)
  lines(x = 1:ncol(trace.m), y = apply(rel.err.m, 1, median), 
        col = "red", lwd = 5)
}
diagnoseBiomass2(lean.fit)
diagnoseBiomass2(obese.fit)








