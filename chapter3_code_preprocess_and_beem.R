library(beemStatic)
library(ggraph)

#MERGE OTU TABLE BY TAXONOMY
library(phyloseq)
library(microbiome)
library(knitr)


path_for_otu <- "C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\otu_table.csv"
path_for_taxa <- "C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\taxonomy.csv"
path_for_meta <- "C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\metadata.csv"

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

#TRANSFORM TO RELATIVE ABUNDANCE
transform_to_rel_ab <- function(phyloseq_object) {
  rel_ab <- transform_sample_counts(phyloseq_object, function(x) x/sum(x))
  return(rel_ab)
}
my_data_rel <- transform_to_rel_ab(my_phylum)
#my_data_rel = transform_sample_counts(my_data, function(x) x/sum(x))

#FILTER ONLY OTUs WITH MEAN GREATER THEN 10E-5
filter_my_data <- function(phyloseq_object) {
  filtered <- filter_taxa(phyloseq_object, function(x) mean(x) > 1e-5, TRUE)
  return(filtered)
}
my_data_fr = filter_taxa(my_data_rel, function(x) mean(x) > 1e-5, TRUE)

# EXTRACT ABUNDANCE MATRIX FROM PHYLOSEQ OBJECT
extract_abundance <- function(phyloseq_object){
  OTU2 = as(otu_table(phyloseq_object), "matrix")
  # transpose if necessary
  if(taxa_are_rows(phyloseq_object)){OTU2 <- t(OTU2)}
  # Coerce to data.frame
  OTUdf2 = as.data.frame(OTU2)
  return(OTUdf2)
}
OTU_phylum_data <- extract_abundance(my_data_fr)

#SEPARATE INTO DBS BASED ON METADATA VALUE (Lean/Obese)
lean_names <- row.names(meta.file)[meta.file$Var == "Lean"]
obese_names <- row.names(meta.file)[meta.file$Var == "Obese"]
View(otu.file)
otu.file.t <- t(otu.file)
mySubset_lean <- otu.file.t[ row.names(otu.file.t) %in% lean_names, ]
mySubset_lean_1 <- mySubset_lean[, colSums(mySubset_lean != 0) > 0]

mySubset_obese <- otu.file.t[ row.names(otu.file.t) %in% obese_names, ]
mySubset_obese_1 <- mySubset_obese[, colSums(mySubset_obese != 0) > 0]

mySubset_lean_meta <- meta.file[row.names(meta.file) %in% lean_names,]
write.csv(mySubset_lean_meta,"C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\otu_table_phylum_lean_meta.csv", row.names = TRUE)

mySubset_obese_meta <- meta.file[row.names(meta.file) %in% obese_names,]
write.csv(mySubset_obese_meta,"C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\otu_table_phylum_obese_meta.csv", row.names = TRUE)


#PREPROCESS SEPARATED DFs
library(funrar)

rel_mat_lean <- make_relative(mySubset_lean_1)
filtered <- filter_taxa(phyloseq_object, function(x) mean(x) > 1e-5, TRUE)



my_phylum <- tax_glom(my_data, taxrank="Phylum")

#PREPARE DATA FOR BEEM NOW
OTU_phylum_data <- t(OTU_phylum_data)
OTU_phylum_data <- OTU_phylum_data[-c(5,6,7),] #Remove low abundance values

df <- OTU_phylum_data

aux.poisson <- function(v, depth){
  v <- round(v/sum(v) * depth)
  sapply(v, function(x) rpois(1, x))
}
#df <- df_pp
df.noise <- apply(df, 2, aux.poisson, depth=0.5)
dat.w.noise=df.noise
dat.wo.noise=df
biomass.true=colSums(df)
beemDemo <- list(dat.w.noise=df.noise, dat.wo.noise=df, biomass.true=colSums(df))

res <- func.EM(dat.w.noise, ncpu=4, scaling=median(biomass.true), max.iter=5)

diagnoseFit(res, dat.w.noise, annotate = TRUE)

diagnoseBiomass(res, true.biomass = biomass.true)

plot(beem2biomass(res), biomass.true, xlab='BEEM biomass estimation', ylab='True biomass')

est <- beem2param(res)

showInteraction(res, dat.w.noise)

#SAVE PHYLUM OTU AS CSV 
write.csv(OTU_phylum_data,"C:\\Users\\Marko\\Desktop\\Ena\\data\\data_to_beem_data\\Turnbaugh\\otu_table_phylum.csv", row.names = TRUE)


#PLOT BY PHYLUM
plot_bar(my_data_fr, fill="Phylum")
