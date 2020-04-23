---
title: "Supplementary Code"
author: "Katariina Pärnänen"
date: "Sep/13/2019"
output:
  github_document: default
---


## Set environment
Acquire required packages. Installation takes less than 15 minutes on a standard laptop computer. Run the code block below to install the required packages. Code run using R version 3.6.1

```{r, echo=TRUE, error=FALSE, message=FALSE}

## ## 
# install.packages("vegan")
## biocLite('phyloseq')
# install.packages("ggplot2")
# install.packages("knitr")
# install.packages("rmarkdown")
# install.packages("cowplot")
# install.packages("viridis")
# install.packages("MASS")
# install.packages("multcomp")
# install.packages("scales")
# install.packages("grid")
# install.packages("reshape2")
# install.packages("randomForest")
# install.packages("dplyr") 
# install.packages("rfUtilities")
# install.packages("Hmisc")
# install.packages("data.table")
# install.packages("ggrepel")
# install.packages("GGally")
# install.packages("jtools")
# install.packages("sjPlot")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# install.packages("lmSupport")
# install.packages("ggsignif")

#Load libraries
## ## 
library(vegan)
library(phyloseq)
library(ggplot2)
library(knitr)
library(DESeq2)
library(rmarkdown)
library(dplyr)
library(cowplot)
library(randomForest)
library(caret)
library(data.table)
library(ggrepel)
library(lme4)
library(lmSupport)
library(svglite)
library(MASS)
library(jtools)
library(multcomp)
library(wesanderson)
library(ggsignif)
##

```

## Read in metaphlan results on species level
The abundance_species.txt file is created from the merged abundance table from Metaphlan. Lines with 's__' are picked and the lines with 't__' are removed. The abundance table then has only species level entries. Column with OTU and running number is added after this.

```{r, warning=FALSE}

metaphlan_sp <- as.matrix(read.table("merged_abundance_table_species.txt", fill= 1, header= T, row.names = 1, check.names = F))
```
## Read in taxonomy table for species level
Taxonomy table is created from the species level merged abundance table with awk and sed scripts.
```{r, warning=FALSE}
tax_sp<- read.table(("tax_table_species.txt"), fill=1, row.names=1)
tax_sp<-apply(tax_sp, 2, function(y) (gsub(".__", "", y)))

```

## Read in metadata
```{r, warning=FALSE}
sample_data<-read.csv(as.matrix("ForAnalysis_NEC_metadata_46.csv"), header=T, row.names = 1,sep = ";", stringsAsFactors=FALSE)
sample_data$Pair<-substring(row.names(sample_data), 1, 6)
sample_data$Sample<-row.names(sample_data)
sample_data$Gest_age2<-cut(sample_data$Gest_age, breaks=3)
sample_data$Antibiotics<-sample_data$Inf_AB=="Y"
sample_data$Donor_Milk<-sample_data$Donor_Breast_Milk=="y"

```

## Merge into a phyloseq object
```{r, warning=FALSE}
PHY_SP <- phyloseq(otu_table(metaphlan_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_sp)), sample_data(sample_data))
```

## Change the taxonomic levels to something meaningful
```{r, warning=FALSE}
colnames(tax_table(PHY_SP)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Remove viruses
PHY_SP<-subset_taxa(PHY_SP, Domain!="Viroids")
PHY_SP<-transform_sample_counts(PHY_SP, function(x) x/sum(x)*100)

# Add dominant genus and family of each sample to the sample data

find.top.taxa <- function(x,taxa){
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- (tax[k,])
  m <- data.frame(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
toptaxa<-find.top.taxa(PHY_SP,"Genus")

toptaxa_fam<-find.top.taxa(PHY_SP,"Family")


top_genus<-toptaxa$Genus
sample_data$Top_genus<-top_genus
sample_data$Top_fam<-toptaxa_fam$Family

# Change all but mentioned to other
sample_data$Top_genus2<-sample_data$Top_genus

sample_data<-sample_data %>%
     mutate(Top_genus2=replace(Top_genus2, !Top_genus%in%c("Escherichia", "Klebsiella", "Veillonella"), "Other")) %>%
     as.data.frame(row.names = row.names(sample_data)) 


# Make phyloseq object again
PHY_SP <- phyloseq(otu_table(metaphlan_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_sp)), sample_data(sample_data))

## Change the taxonomic levels to something meaningful
colnames(tax_table(PHY_SP)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Remove viruses
PHY_SP<-subset_taxa(PHY_SP, Domain!="Viroids")
# Transform sample counts back to 100 %
PHY_SP<-transform_sample_counts(PHY_SP, function(x) x/sum(x)*100)

```
## Read in Metaxa2 results

```{r, results="hide", messages =FALSE, warnings =FALSE}
# Read in metaxa OTU table and tax table
metaxa <- as.matrix(read.table("metaxa_mat.txt", fill= 1, header= T, row.names = 1, check.names = F))
metaxa_tax <- read.table(("tax_table_metaxa.txt"), fill=1, row.names=1, header=F)
colnames(metaxa)<-gsub("_", "-", colnames(metaxa))

# Make phyloseq object
PHY_mtx <- phyloseq(otu_table(metaxa, taxa_are_rows = TRUE), tax_table(as.matrix(metaxa_tax)), sample_data(sample_data))

# Change the levels to something meaningful
colnames(tax_table(PHY_mtx)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Exclude Eukaryota and unknown domains and samples with less than 20 16S SSU matches
PHY_mtx_mod <- subset_taxa(PHY_mtx, !Domain%in%c("Eukaryota", "unknown"))
PHY_mtx<- subset_samples(PHY_mtx_mod, sample_sums(PHY_mtx)>20)

# Turn into relative data
PHY_mtx_rel<-transform_sample_counts(PHY_mtx, function(x) x/sum(x))
toptaxa_mtx<-find.top.taxa(PHY_mtx_rel, "Family")
sample_data$Top_genus_mtx<-toptaxa_mtx$Family

sample_data(PHY_SP)<-sample_data

```

# Analysis of ARGs and MGEs
## Read in ARG and MGE bowtie2 mapping results

```{r, warning=FALSE}
ARG_bt<-as.matrix(read.table("ARG_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths <- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_bt_res<-as.matrix(read.table("ARG_genemat_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))

MGE_bt<-as.matrix((read.table("MGE_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F)))
MGE_tax<-read.table("/Users/kparnane/Documents/Maito/Nextseq_milk1-2/MGE_tax_table.txt", fill=1, row.names = 1, header=F)


# Divide by genelengths
# Divide by ARG gene lengths
ARG_bt<-ARG_bt[row.names(ARG_bt)%in%row.names(ARG_lengths),]
ARG_bt_res<-ARG_bt_res[row.names(ARG_bt_res)%in%row.names(ARG_lengths_res),]
arg_length_norm <- ARG_bt/ARG_lengths[,1]
arg_length_norm_res <- ARG_bt_res/ARG_lengths_res[,1]


# Read in SSU_counts, this is always > 1000
SSU_counts<-read.table(as.matrix("SSU_counts.txt"), header=F, row.names = 1,sep = " ", stringsAsFactors=FALSE)
SSU_counts[,1]<-as.numeric(SSU_counts[,1])
ARG_SSU_length_norm<-t(t(arg_length_norm)/SSU_counts[,1])*1541
ARG_SSU_length_norm_res<-t(t(arg_length_norm_res)/SSU_counts[,1])*1541

# Check if the division is correct with the first column
identical(ARG_SSU_length_norm[,1], arg_length_norm[,1]/22955*1541)

# Repeat for MGE sequences

MGE_length <- as.matrix(read.table("/Users/kparnane/Documents/Maito/Nextseq_milk1-2/MGE_genelenghts.txt", fill= 1, header= F, row.names = 1, check.names = F))

MGE_bt<-MGE_bt[row.names(MGE_bt)%in%row.names(MGE_length),]
mge_length_norm <- MGE_bt/MGE_length[,1]
MGE_SSU_length_norm<-t(t(mge_length_norm) / SSU_counts[,1])*1541
MGE_SSU_length_norm[is.na(MGE_SSU_length_norm)]<-0

# Remove too long and too short genes
temp<-MGE_length[(MGE_length[,1] < 4000) & (MGE_length[,1] > 100), ]
MGE_SSU_length_norm<-MGE_SSU_length_norm[row.names(MGE_SSU_length_norm)%in%c(names(temp)),]


```

## Make phyloseq object

```{r, results="hide", messages =FALSE, warnings =FALSE}
# Make phyloseq objects from ARG and MGE matrixes

ARG_res<-phyloseq(otu_table(ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(sample_data), tax_table(as.matrix(ARG_tax)))

MGE<-phyloseq(otu_table(MGE_SSU_length_norm, taxa_are_rows = T), sample_data(sample_data), tax_table(as.matrix(MGE_tax)))
MGE_tax$V2<-gsub('_', ' ', MGE_tax$V2)
MGE_tax$V2<-gsub('insertion element ', '', MGE_tax$V2)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)

# Remove the genes which have no matches and genes annotated as Bifidobacteria
MGE_mod<-otu_table(MGE)[rowSums(otu_table(MGE))>0]
MGE_PHY <- phyloseq(MGE_mod, sample_data(sample_data), tax_table(as.matrix(MGE_tax)))

# Using only resfinder hits remove genes which have no matches and make phyloseq object
ARG_mod2_res<-otu_table(ARG_res)[rowSums(otu_table(ARG_res))>0]
ARG_PHY_res<-phyloseq(ARG_mod2_res, sample_data(sample_data), tax_table(as.matrix(ARG_tax)))

# Remove bad gene
ARG_PHY_res2<-(subset_taxa(ARG_PHY_res, !V3%in%c("oqxB")))

# Make object also for not normalized data
ARG_PHY_unnorm<-phyloseq(otu_table(ARG_bt_res, taxa_are_rows = TRUE), sample_data(sample_data), tax_table(as.matrix(ARG_tax)))
ARG_PHY_unnorm<-subset_taxa(ARG_PHY_unnorm, !V3%in%c("oqxB"))

# Add 16S counts
sample_data(ARG_PHY_unnorm)$SSU_counts<-SSU_counts[,1]
```

# Humann2

```{r, results="hide", messages =FALSE, warnings =FALSE}
human_sp <- as.matrix(read.table("all_genefamilies_cpmEC_genemat.tsv", fill= 1, header= T, check.names = F, sep = "\t"))
OTU<- paste("OTU", 1:nrow(human_sp), sep="")
rownames(human_sp)<- OTU
```

## Read in taxonomy table of genes

```{r, warning=FALSE}
tax_EC<- read.table(("EC_taxtable.txt"), fill=TRUE, header=TRUE, sep ="\t", row.names = 1, quote =  "")
OTU<- paste("OTU", 1:nrow(tax_EC ), sep="")
rownames(tax_EC )<- OTU

```

## Merge into a phyloseq object

```{r, warning=FALSE}
PHY_humann <- phyloseq(otu_table(human_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_EC)), sample_data(sample_data))
```


## Ordinate and calculate diversities
```{r, warning=FALSE} 


## ## ## ## ## ## ## Pairwise adonis function## ## ## ## ## 

## start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'horn', p.adjust.m ='fdr')
{
library(vegan)
  

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()


for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
return(pairw.res)

} 

## end copy here

## ## ## ## ## ## ## ## ## # 

# Metaphlan

pairwise.adonis(x=sqrt(t(otu_table(PHY_SP))), factors=sample_data(PHY_SP)$Top_genus2, sim.function='vegdist', sim.method='horn',p.adjust.m='fdr')

PHY_SP_sqrt<-phyloseq(sqrt(otu_table(PHY_SP, taxa_are_rows = TRUE)), sample_data(sample_data(PHY_SP)), tax_table(PHY_SP))
# Ordinate Metaphlan results
PHY_SP_ord<-ordinate(PHY_SP_sqrt, method="PCoA", distance="horn")

mylabels <- c(expression(paste(italic("Escherichia"))), 
              expression(paste(italic("Klebsiella"))), 
               expression(paste(italic("Veillonella"))), 
               "Other") 

p<-plot_ordination(PHY_SP_sqrt, PHY_SP_ord, color="Top_genus2") + guides(color=guide_legend(ncol=2))
p$layers <- p$layers[-1]
metaphlan.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", breaks=c("Escherichia", "Klebsiella", "Veillonella", "Other"), labels=mylabels) +theme_classic() + geom_point(size=4, alpha=0.7)+ labs(title="Microbiota") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) + stat_ellipse(level = 0.90, linetype=2, show.legend = FALSE) +
  geom_point(shape=21, size=4, color="black") +
  ylab("Axis 2 (19.4 %)") +
  xlab("Axis 1 (25.9 %)")+
 # ylim(-0.75, 0.75) +
 # xlim(-0.75, 0.75) +
  coord_equal() 

# Check are the separations significant for ARGs

pairwise.adonis(x=sqrt(t(otu_table(ARG_PHY_res2))), factors=sample_data(ARG_PHY_res2)$Top_genus2, sim.function='vegdist', sim.method='horn',p.adjust.m='fdr')

# Check separation in MGEs

pairwise.adonis(x=sqrt(t(otu_table(MGE_PHY))), factors=sample_data(MGE_PHY)$Top_genus2, sim.function='vegdist', sim.method='horn',p.adjust.m='fdr')

# Ordinate ARG results

PHY_ARG_sqrt<-phyloseq(sqrt(otu_table(ARG_PHY_res2, taxa_are_rows = TRUE)), sample_data(sample_data(ARG_PHY_res2)), tax_table(ARG_PHY_res2))

PHY_ARG_ord<-ordinate(PHY_ARG_sqrt, method="PCoA", distance="horn")

# Plot

p<-plot_ordination(PHY_ARG_sqrt, PHY_ARG_ord, color="Top_genus2")
p$layers <- p$layers[-1]
ARG.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", breaks=c("Escherichia", "Klebsiella", "Veillonella", "Other")) + geom_point(size=4, alpha=0.7)  +theme_classic() + labs(title="ARGs") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, linetype=2, show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
  ylab("Axis 2 (11.8 %)") +
  xlab("Axis 1 (23.1 %)")
  # ylim(-0.85, 0.85) +
  # xlim(-0.85, 0.85) +

# Ordinate MGE results

PHY_MGE_sqrt<-phyloseq(sqrt(otu_table(MGE_PHY, taxa_are_rows = TRUE)), sample_data(sample_data(MGE_PHY)), tax_table(MGE_PHY))

PHY_MGE_ord<-ordinate(PHY_MGE_sqrt, method="PCoA", distance="horn")

# Plot

p<-plot_ordination(PHY_MGE_sqrt, PHY_MGE_ord, color="Top_genus2")
p$layers <- p$layers[-1]
MGE.plot<-p+scale_color_brewer(palette = "BrBG", "Dominant genus", breaks=c("Escherichia", "Klebsiella", "Veillonella", "Other")) + geom_point(size=4, alpha=0.7)   +theme_classic() + labs(title="MGEs") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, linetype=2, show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
  ylab("Axis 2 (14.8 %)") +
  xlab("Axis 1 (18.2 %)")
  # ylim(-0.85, 0.85) +
  # xlim(-0.85, 0.85) +


# Ordinate Humann2 results

PHY_humann <- phyloseq(otu_table(human_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_EC)), sample_data(sample_data))
PHY_humann_gene_prop<-transform_sample_counts(PHY_humann, function(x) x/sum(x))
PHY_humann_gene_prop_sqrt <- phyloseq(sqrt(otu_table(PHY_humann_gene_prop, taxa_are_rows = TRUE)), sample_data(sample_data(PHY_humann_gene_prop)), tax_table(PHY_humann_gene_prop))
PHY_humann_ord<-ordinate(PHY_humann_gene_prop_sqrt, method="PCoA", distance="horn")


# Plot
p<-plot_ordination(PHY_humann_gene_prop_sqrt, PHY_humann_ord, color="Top_genus2")
p$layers <- p$layers[-1]
humann.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", breaks=c("Escherichia", "Klebsiella", "Veillonella", "Other")) + geom_point(size=4, alpha=0.7)  +
  theme_classic() + labs(title="Functional genes") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, linetype=2, legend=FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
  ylab("Axis 2 (20.1 %)")+
  xlab("Axis 1 (32.4 %)")
   # ylim(-0.5, 0.5) +
  # xlim(-0.5, 0.5) +
  coord_equal()


# Check significance of Humann clusters
pairwise.adonis(sqrt(t(otu_table(PHY_humann_gene_prop))), factors = sample_data(PHY_humann_gene_prop)$Top_genus2, sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")


# Diversity
df <- data.frame(DIV=diversity(t(otu_table(PHY_SP)), index = "shannon"),
                 DIV2=colSums(otu_table(PHY_SP)>0),
                  Top_genus2=sample_data(PHY_SP)$Top_genus2, Top_genus2=sample_data(PHY_SP)$Top_genus2)

metaphlan.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(alpha=FALSE)+
theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("Microbiota")

# Shannon diversity 
a0<-aov(diversity(t(otu_table(PHY_SP)), index = "shannon")~sample_data(PHY_SP)$Top_genus2, data=as.data.frame(otu_table(PHY_SP)))
TukeyHSD(a0)


df <- data.frame(DIV=diversity(t(otu_table(PHY_humann_gene_prop)), index = "shannon"), DIV2=colSums(otu_table(PHY_humann_gene_prop)>0), Top_genus2=sample_data(PHY_humann_gene_prop)$Top_genus2, Top_genus2=sample_data(PHY_humann_gene_prop)$Top_genus2)
humann.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("Functional genes")

a0<-aov(diversity(t(otu_table(PHY_humann_gene_prop)), index = "shannon")~sample_data(PHY_humann_gene_prop)$Top_genus2, data=as.data.frame(otu_table(PHY_humann_gene_prop)))
TukeyHSD(a0)

# Calculate diverisity and make plot

df <- data.frame(DIV=diversity(t(otu_table(ARG_PHY_res2)), index = "shannon"),
                 DIV2=colSums(otu_table(ARG_PHY_res2)>0),
                 Top_genus2=sample_data(ARG_PHY_res2)$Top_genus2, Top_genus2=sample_data(ARG_PHY_res2)$Top_genus2)

ARG.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank())  +labs(y="Shannon diversity", x='') + ggtitle("ARGs")

# Analysis of variance

a0<-aov(diversity((t(otu_table(ARG_PHY_res2))), index = "shannon")~sample_data(ARG_PHY_res2)$Top_genus2, data=as.data.frame(otu_table(ARG_PHY_res2)))
TukeyHSD(a0)

# Calculate diverisity and make plot

df <- data.frame(DIV=diversity(t(otu_table(MGE_PHY)), index = "shannon"),
                 DIV2=colSums(otu_table(MGE_PHY)>0),
                 Top_genus2=sample_data(MGE_PHY)$Top_genus2, Top_genus2=sample_data(MGE_PHY)$Top_genus2)

MGE.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG","Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("MGEs")

# Analysis of variance
a0<-aov(diversity((t(otu_table(MGE_PHY))), index = "shannon")~sample_data(MGE_PHY)$Top_genus2, data=as.data.frame(otu_table(MGE_PHY)))
TukeyHSD(a0)

# Cowplot

legend<-get_legend(plot=metaphlan.plot)

cowplot::plot_grid(metaphlan.plot +theme(legend.position = "none"),
                           humann.plot+theme(legend.position = "none"),
                           ARG.plot+theme(legend.position = "none"), 
                           MGE.plot+theme(legend.position = "none"),
                           labels = c("auto"), align = "hv", ncol = 2)
```

## Random forests for data exploration to pick best predictors for sum abundances
```{r, results="hide", messages =FALSE, warnings =FALSE}
# Save mean 16S read count in samples for normalization
mean16S_counts<-mean(SSU_counts[,1])

# Random forests to get the best predictors for ARG sum abundance
# The dataframe has been simplified from the whole set of metadata to demonstrate the use of RFs
df <- data.frame(SUM=((sample_sums(otu_table(ARG_PHY_res)))),
                 MGE_SUM=((sample_sums(otu_table(MGE_PHY)))),
                 Delivery_met=sample_data(ARG_PHY_res)$Delivery.Method,
                Diet=sample_data(ARG_PHY_res2)$DIET_CODES_4,
                 Inf_AB=sample_data(ARG_PHY_res)$Inf_AB,
                 Mat_AB=sample_data(ARG_PHY_res)$Mat_AB,
                 M.Seqs=sample_data(ARG_PHY_res)$M.Seqs,
                 Gest_age=sample_data(ARG_PHY_res)$Gest_age,
                 Age_sampling=sample_data(ARG_PHY_res2)$Age_sampling,
                 Inf_inf=sample_data(ARG_PHY_res)$Inf_Inf)


# Save variables for ARG model
y<-df$SUM
x<-df[,3:10]
dfr<-cbind(x, y)

# Train models
fitControl <- trainControl(method="cv")
rfModel <- train(x=x,y=y,data=data, method="rf", rControl=fitControl,importance=TRUE)

# Get importance of variables
varImp(rfModel)


```
## Build model for ARG and MGE sum abundance using gamma distributed GLMs

```{r, results="hide", messages =FALSE, warnings =FALSE}
# ARGs
# Get integers for negative binomial distribution by multiplying the relative gene counts to 16S with mean 16S counts 

mean16S_counts<-mean(SSU_counts[,1])

#Add the explanatory variables to a dataframe
df <- data.frame(MGE_SUM=(mean16S_counts*(sample_sums(otu_table(MGE_PHY)))),
                SUM=((sample_sums(otu_table(ARG_PHY_unnorm)))), 
                SUM2=sample_sums(ARG_PHY_res),
                ENTERO=sample_sums(subset_taxa(PHY_SP, Family=="Enterobacteriaceae")), 
                GAMMA=sample_sums(subset_taxa(PHY_SP, Class=="Gammaproteobacteria")),
                SSU_counts=SSU_counts[,1],
                 Gest_Age=sample_data(ARG_PHY_res)$Gest_age,
                  DIET_CODES_4=sample_data(ARG_PHY_res)$DIET_CODES_4,
                DIET_CODES_1=sample_data(ARG_PHY_res)$DIET_CODES_1,
                BOV=sample_data(ARG_PHY_res)$DIET_bov,
                 Inf_Inf=sample_data(ARG_PHY_res)$Inf_Inf,
                M.Seqs=sample_data(ARG_PHY_res)$M.Seqs,
                Twin=sample_data(ARG_PHY_res)$Twin,
                 Inf_AB=sample_data(ARG_PHY_res)$Inf_AB,
                Suspected_GI_problem=sample_data(MGE_PHY)$Suspected_GI_problem,
                 Gest_age=sample_data(ARG_PHY_res)$Gest_age,
                Mat_AB=sample_data(MGE_PHY)$Mat_AB,
                Age_sampling=sample_data(ARG_PHY_res2)$Age_sampling,
                Formula=sample_data(MGE_PHY)$Formula,
                Fortifier=sample_data(MGE_PHY)$Milk_Fortifier,
                Donor=sample_data(MGE_PHY)$Donor_Breast_Milk,
                Mom=sample_data(MGE_PHY)$Mom_Breast_Milk,
                Formula_perc=sample_data(ARG_PHY_res2)$Week_1.Formula+sample_data(ARG_PHY_res2)$Week_2.Formula,
                Fortifier_perc=sample_data(ARG_PHY_res2)$Week_1..Fortified_breast_milk+sample_data(ARG_PHY_res2)$Week_2_Fortified_breast_milk,
                FF_perc=sample_data(ARG_PHY_res2)$Week_1.Formula+sample_data(ARG_PHY_res2)$Week_2.Formula+sample_data(ARG_PHY_res2)$Week_1..Fortified_breast_milk+sample_data(ARG_PHY_res2)$Week_2_Fortified_breast_milk,
                Milk_perc=sample_data(ARG_PHY_res2)$Week_1.Breast_milk+sample_data(ARG_PHY_res2)$Week_2.Breast_milk,            Pair=sample_data(ARG_PHY_res)$Pair,
                week1BM=sample_data(ARG_PHY_res2)$Week_1.Breast_milk,
                week2BM=sample_data(ARG_PHY_res2)$Week_2.Breast_milk,
                Delivery_mode=sample_data(ARG_PHY_res2)$Delivery.Method,
                Top_fam=sample_data(ARG_PHY_res2)$Top_fam)



#Check also with older infants, since formula was mostly given to older infants.
df2<-df[df$Gest_Age>32,]

#Check also with younger infants, since formula was mostly given to older infants.
df3<-df[df$Gest_Age<32,]

#Check also for infants who have formula
df4<-df[df$Formula=="y",]
df4_2<-df[df$Formula=="n",]

df$AB_exp=(df$Inf_AB=="Y")+(df$Mat_AB=="Y")
df$Milk_exp=(df$week1BM)+(df$week2BM)>100


#Check also by removing the two super high values


##################################

#Test if there is overdispersion

mp <- glm(SUM~Gest_age+DIET_CODES_4, family=poisson, data=df)

summary(mp)



qchisq(0.95, df.residual(mp))

deviance(mp)

pr <- residuals(mp,"pearson")

sum(pr^2)

#95% cutoff value for chisq
phi <- sum(pr^2)/df.residual(mp)


round(c(phi,sqrt(phi)),4)


# Quasipoisson

mq <- glm(SUM~Gest_age+DIET_CODES_4, family=quasipoisson, data=df)

summary(mq)
se <- function(model) sqrt(diag(vcov(model)))

round(data.frame(p=coef(mp), q=coef(mq), se.p=se(mp), se.q=se(mq), ratio=se(mq)/se(mp)), 4)

# Negative binomial

mnb <- glm.nb(SUM~Gest_age+DIET_CODES_4, data=df)

summary(mnb)

# Variance
1/mnb$theta

mnbg <- glm(SUM~Gest_age+DIET_CODES_4, family=negative.binomial(mnb$theta), data=df)

all(abs(coef(mnb)==coef(mnbg)) < 1e-12)


v = 1/mnb$theta
qgamma((1:3)/4, shape = 1/v, scale = v)

round(data.frame(p=coef(mp),q=coef(mq),nb=coef(mnb), se.p=se(mp),se.q=se(mq),se.nb=se(mnb)),4)

deviance(mnbg)


######################################
#Build the models for ARG abundance with overdispersion

# Gamma distrubuted GLM with log link
fit1_gamma<-glm(SUM2~DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit1_gamma)
summ(fit1_gamma, exp=TRUE)


# Tukey
glht.mod <- glht(fit1_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))


fit2_gamma<-glm(SUM2~Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit2_gamma)
summ(fit2_gamma, exp=TRUE)

fit2_gamma_MGE<-glm(MGE_SUM~Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit2_gamma)
summ(fit2_gamma, exp=TRUE)

glht.mod <- glht(fit2_gamma_MGE, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))


qchisq(0.95, df.residual(fit2_gamma))

deviance(fit2_gamma)

cor(model.matrix(fit2_gamma)[,-1])


glht.mod <- glht(fit2_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))


fit2_gamma_ntw<-glm(SUM2~Gest_age+DIET_CODES_4, data=df[df$Twin=="n",], family=Gamma(link="log"))
summary(fit2_gamma_ntw)
summ(fit2_gamma_ntw)

fit2_gamma_nof<-glm(SUM2~Gest_age+DIET_CODES_4, data=df[df$Formula=="n",], family=Gamma(link="log"))
summary(fit2_gamma_nof)
summ(fit2_gamma_nof)



library(pscl)

pchisq(summary(fit2_gamma)$deviance, 
           summary(fit2_gamma)$df.residual
           )


fit3_gamma<-glm(SUM2~Gest_age*DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit3_gamma)
summ(fit3_gamma)
pR2(fit3_gamma)

temp<-cbind(df, 
      Mean = predict(fit3_gamma, newdata=df, type="response"), 
      SE = predict(fit3_gamma, newdata=df, type="response", se.fit=T)$se.fit
      )

fit4_gamma<-glm(SUM2~Inf_AB+Gest_age*DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit4_gamma)
summ(fit4_gamma, exp = TRUE)

fit5_gamma <- glm(SUM2~Age_sampling+Inf_AB+Gest_age*DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit5_gamma)

fit6_gamma <- glm(SUM2~Age_sampling+Mat_AB+Gest_age*DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit6_gamma)

# Perform chisquared test

anova_NEC<-anova(fit1_gamma, fit2_gamma, fit3_gamma, fit4_gamma, fit5_gamma, test="Chisq")


# Use Tukey's post hoc test for the model 2
glht.mod <- glht(fit2_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))


fit5_gamma<-glm(SUM2~Gest_age+DIET_CODES_1, data=df, family=Gamma(link="log"))
summary(fit5_gamma)

glht.mod <- glht(fit5_gamma, mcp(DIET_CODES_1="Tukey"))
summary(glht(glht.mod))


library(ggsignif)

annotation_df <- data_frame(start=c("FF", "EBM", "EBM"), 
                            end=c("F", "F", "FF"),
                            y=c(0.5, 0.75, 1),
                            label=c("**", "**", "ns"))

# Save ARG.sum.plot

ARG.sum.plot<-ggplot(df, aes(x=DIET_CODES_4, y=log10(SUM2))) + 
    scale_fill_brewer(palette="BrBG")+
    geom_line() +
    geom_boxplot()+
    geom_jitter(data=df, aes(x=DIET_CODES_4, y=log10(SUM2), color=DIET_CODES_4, fill=DIET_CODES_4), size=5, width=0.3, shape=21, color="black")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, vjust=0.65, hjust=0.8), text=element_text(size=16)) +
    labs(y="Log10 sum abundance/16S rRNA gene", x='') +
    guides(fill=FALSE, alpha=FALSE)+
    labs(title="ARGs") +
   geom_signif(data=annotation_df,
              aes(xmin=start, xmax=end, annotations=label, y_position=y),
              textsize = 4, vjust = -0.2,
              manual=TRUE) +
    scale_x_discrete(breaks=c("EBM","F","FF"),
        labels=c("Breast milk", "Formula", "Fortifier"))

ARG.sum.plot


# Plot effect of gestational age

gamma_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "Gamma"(link="log")))}


gest.age.plot<-ggplot(df, aes(x=Gest_age, y=(SUM2), fill=Formula, color=Formula)) +
  gamma_smooth()+
  geom_point(shape=21, color="black", size=5)+
    theme_classic()+
    scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes") ) +
     scale_color_manual(values=c("darkgoldenrod3", "azure4"),labels=c("No", "Yes"))+
    labs(x="Gestational age (weeks)", y="ARG sum abundance/\n16S rRNA gene") +
  theme(axis.title=element_text(size=16), axis.text.x = element_text(size=rel(1.2)), plot.title = element_text(size=16), axis.text=element_text(size=10))+
  scale_y_sqrt()+
       guides(fill = guide_legend(override.aes = list(size = 2)),
    size = guide_legend(override.aes = list(linetype = 0)))


gest.age.plot

# Build  model for MGE abundace

fit1<-glm((MGE_SUM)/mean16S_counts~Formula, data=df, family="Gamma"(link="log"))
summary(fit1)

# Plot
annotation_df <- data.frame(start=c("FF", "EBM", "EBM"), 
                            end=c("F", "F", "FF"),
                            y=c(1.5, 1.75, 2),
                            label=c("*", ".", "ns"))

MGE.sum.plot<-ggplot(df, aes(x=DIET_CODES_4, y=log10(MGE_SUM/mean16S_counts))) + 
    scale_fill_brewer(palette = "BrBG") +
    geom_line() +
  geom_boxplot()+
   geom_jitter(data=df, aes(x=DIET_CODES_4, y=log10(MGE_SUM/mean16S_counts), color=DIET_CODES_4, fill=DIET_CODES_4), size=5, width=0.3, shape=21, color="black")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45, vjust=0.65, hjust=0.8), text=element_text(size=16)) +
   geom_signif(data=annotation_df,
              aes(xmin=start, xmax=end, annotations=label, y_position=y),
              textsize = 4, vjust = -0.2,
              manual=TRUE) +
    labs(y="Log10 sum abundance/16S rRNA gene", x='') +
    guides(fill=FALSE, alpha=FALSE)+
    labs(title="MGEs") +
   scale_x_discrete(breaks=c("EBM","F","FF"),
        labels=c("Breast milk", "Formula", "Fortifier"))

MGE.sum.plot


# Build  model for Enterobacteriaceae abundace
fit1<-glm(ENTERO/100~Gest_age+Formula, data=df, family="quasibinomial")
summary(fit1)


# Build model for Gammaproteobacteria abundance
fit1_gamma<-glm(GAMMA/100~Gest_age+Formula, data=df, family="quasibinomial")
summary(fit1_gamma)


# Make cowplot
cowplot:::plot_grid(ARG.sum.plot, MGE.sum.plot, labels="auto", align="h", nrow=1)


```


# DESEQ

DESEQ analysis for ARGs 

```{r, warning=FALSE, message=FALSE}
# ARGs and breastfeeding exclusively or FF vs formula
# Pseudocounts
temp_arg<-ARG_SSU_length_norm_res[,]*5*10^4+1
# Change the gene names to be harmonious in synthax
ARG_tax2<-ARG_tax
ARG_tax2$V3<-gsub('Erm', 'erm', ARG_tax2$V3)
ARG_tax2$V3<-gsub('erm\\(C\\)', 'ermC', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ErmC', 'ermC', ARG_tax2$V3)
ARG_tax2$V3<-gsub('erm\\(B\\)', 'ermB', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ErmB', 'ermB', ARG_tax2$V3)
ARG_tax2$V3<-gsub('Cfx', 'cfx', ARG_tax2$V3)
ARG_tax2$V3<-gsub('APH', 'aph', ARG_tax2$V3)
ARG_tax2$V3<-gsub('FosA', 'fosA', ARG_tax2$V3)
ARG_tax2$V3<-gsub('AAC', 'aac', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ANT', 'ant', ARG_tax2$V3)


# Make phyloseq object for DESEQ2
ARG_DSQ<-phyloseq(otu_table(temp_arg, taxa_are_rows=T), sample_data(sample_data), tax_table(as.matrix(ARG_tax2)))
ARG_DSQ_glom<-tax_glom(ARG_DSQ, taxrank="V3")
dds_arg = phyloseq_to_deseq2(ARG_DSQ_glom, ~ Formula)
dds_arg$Formula <- relevel(dds_arg$Formula, "n")
dds_arg = DESeq(dds_arg, fitType = "mean", test ="Wald", betaPrior = TRUE)
res_arg = results(dds_arg, cooksCutoff = FALSE)
alpha = 0.05
sigtab_arg = res_arg[which(res_arg$padj < alpha), ]
sigtab_arg = cbind(as(sigtab_arg, "data.frame"), as(tax_table(ARG_DSQ)[rownames(sigtab_arg), ], "matrix"))

# Save data on how many samples the ARG was found in
otu_table(ARG_DSQ)[otu_table(ARG_DSQ)==1]<-0
otu_table(ARG_DSQ)[otu_table(ARG_DSQ)>0]<-1
n<-rowSums(otu_table(ARG_DSQ))

# Make result table
sigtab_arg=merge(sigtab_arg, as.data.frame(n), by=0)

# Filter by prevalence
sigtab_arg<-sigtab_arg[sigtab_arg$n>6,]

# Filter by abundance
sigtab_arg<-sigtab_arg[sigtab_arg$baseMean>4,]

# Save Rrsistance class order
x = tapply(sigtab_arg$log2FoldChange, sigtab_arg$V2, function(x) max(x))
x = sort(x, TRUE)
sigtab_arg$V2 = factor(as.character(sigtab_arg$V2), levels=names(x))
# Gene order
x = tapply(sigtab_arg$log2FoldChange, sigtab_arg$V3, function(x) max(x))
x = sort(x, TRUE)
sigtab_arg$V3 = factor(as.character(sigtab_arg$V3), levels=names(x))

# Plot
a<-ggplot(sigtab_arg, aes(x=V3, y=log2FoldChange, fill=V2, size=baseMean)) + geom_point(alpha=0.7, shape=21, color="black")  +theme_minimal()+theme(axis.text.x = element_text(angle =45, hjust = 0.9, vjust = 0.97, size=rel(0.8), face = "italic"), axis.text.y=element_text(size=rel(0.9)), axis.title.x = element_blank(), axis.title.y = element_text(size=rel(0.8)), legend.text= element_text(size=rel(0.8)), legend.title = element_text(size=rel(0.6)), legend.key.width=unit(0.1, "lines"), legend.key.height =unit(0.5, "lines")) + scale_fill_brewer(palette="BrBG") + scale_size(range=c(3, 20),breaks=c(3, 5, 10, 15, 20)) + labs(fill="") +guides(size=FALSE) + ylim(-5, 7.5) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") + ggtitle("ARGs") +
  guides(fill = guide_legend(override.aes = list(size = 2)))

# MGEs and breastfeeding exclusively vs formula
# Pseudocounts
temp_arg<-MGE_SSU_length_norm[,]*5*10^4+1
MGE_tax2<-MGE_tax
# Change the gene names to be harmonious in synthax
MGE_tax2$V3<-gsub('delta', 'Δ', MGE_tax2$V3)
MGE_tax2$V2<-gsub('delta', 'Δ', MGE_tax2$V2)
MGE_tax2$V2<-gsub('tniA', 'transposase', MGE_tax2$V2)
MGE_tax2$V2<-gsub('tniB', 'transposase', MGE_tax2$V2)
MGE_tax2$V2<-gsub('IS[[:digit:]]+', 'insertion element', MGE_tax2$V2)
MGE_tax2$V2<-gsub('istB', 'IS21 transposase', MGE_tax2$V2)
MGE_tax2$V2<-gsub('istA[[:digit:]]+', 'IS21 transposase', MGE_tax2$V2)
MGE_tax2$V2<-gsub('ISSfl3+', 'insertion element', MGE_tax2$V2)

MGE_DSQ<-phyloseq(otu_table(temp_arg, taxa_are_rows=T), sample_data(sample_data), tax_table(as.matrix(MGE_tax2)))
dds_mge = phyloseq_to_deseq2(MGE_DSQ, ~ Formula)
dds_mge$DIET_CODES_4 <- relevel(dds_mge$Formula, "n")
dds_mge = DESeq(dds_mge, fitType = "mean", test ="Wald", betaPrior = FALSE)
res_mge = results(dds_mge, cooksCutoff = FALSE)
alpha = 0.05
sigtab_mge = res_mge[which(res_mge$padj < alpha), ]
sigtab_mge = cbind(as(sigtab_mge, "data.frame"), as(tax_table(MGE_DSQ)[rownames(sigtab_mge), ], "matrix"))


otu_table(MGE_DSQ)[otu_table(MGE_DSQ)==1]<-0
otu_table(MGE_DSQ)[otu_table(MGE_DSQ)>0]<-1
n<-rowSums(otu_table(MGE_DSQ))

sigtab_mge=merge(sigtab_mge, as.data.frame(n), by=0)

sigtab_mge<-sigtab_mge[sigtab_mge$n>6,]

# Make a figure with ggplot
# MGE class order
x = tapply(sigtab_mge$log2FoldChange, sigtab_mge$V2, function(x) max(x))
x = sort(x, TRUE)
sigtab_mge$V2 = factor(as.character(sigtab_mge$V2), levels=names(x))
# Gene order
x = tapply(sigtab_mge$log2FoldChange, sigtab_mge$V3, function(x) max(x))
x = sort(x, TRUE)
sigtab_mge$V3 = factor(as.character(sigtab_mge$V3), levels=names(x))

# Plot

pal <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000", "#446455", "#FDD262")

b<-ggplot(sigtab_mge, aes(x=V3, y=log2FoldChange, fill=V2, size=baseMean)) + geom_point(alpha=0.9, shape=21, color="black")  +theme_minimal()+theme(axis.text.x = element_text(angle =45, hjust = 0.9, vjust = 0.97, size=rel(0.8), face = "italic"), axis.text.y=element_text(size=rel(0.9) ), axis.title.x = element_blank(), axis.title.y = element_text(size=rel(0.8)), legend.text= element_text(size=rel(0.8)), legend.title = element_text(size=rel(0.6)), legend.key.width=unit(0.1, "lines"), legend.key.height =unit(0.5, "lines")) + scale_fill_manual(values=pal) + scale_size(range=c(3, 20),breaks=c(3, 5, 10, 15, 20)) + labs(fill="") +guides(size=FALSE) + ylim(-10, 10) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") + ggtitle("MGEs") +
  guides(fill = guide_legend(override.aes = list(size = 2)))

# Cowplot

cowplot::plot_grid(a,b, ncol = 1, align ="v")
```

## Mantel's test for correlations between taxa, ARGs and MGEs

```{r, results="hide", messages =FALSE, warnings =FALSE}
#Do a mantel test for ARG and MGE distance matrixes

MGE_dist<-vegdist(t(as.matrix(otu_table(MGE_PHY))), method="horn")
ARG_dist<-vegdist(t(as.matrix(otu_table(ARG_PHY_res2))), method="horn")
mantel(ARG_dist, MGE_dist, method="kendall")

PHY_dist<-vegdist(t(as.matrix(otu_table(PHY_SP))), method="horn")
ARG_dist<-vegdist(t(as.matrix(otu_table(ARG_PHY_res2))), method="horn")
mantel(ARG_dist, PHY_dist, method="kendall")

PHY_dist<-vegdist(t(as.matrix(otu_table(PHY_SP))), method="horn")
MGE_dist<-vegdist(t(as.matrix(otu_table(MGE_PHY))), method="horn")
mantel(MGE_dist, PHY_dist, method="kendall")
```

# Meta-analysis datasets

The following code blocks for individual publications are for compiling the metagenome analysis data from ARG mapping and meta-data to a suitable format for meta-analysis GLM analyses and other analyses.

# Bäckhed
https://doi.org/10.1016%2Fj.chom.2015.04.004
```{r, echo=FALSE, error=FALSE, message=FALSE}

Backhed_ARG_bt_res<-as.matrix(read.table("Backhed_ARG_genemat_resfinder_new.txt", fill= 1, header= T, row.names = 1, check.names = F))


ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Backhed_ARG_bt_res<-Backhed_ARG_bt_res[row.names(Backhed_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Backhed_arg_length_norm_res <- Backhed_ARG_bt_res/ARG_lengths_res[,1]

Backhed_ARG_bt<-Backhed_ARG_bt_res[row.names(Backhed_ARG_bt_res)%in%row.names(ARG_lengths),]

#Read in SSU_counts
Backhed_SSU_counts <- as.matrix(read.table("Backhed_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))



#Read in SSU counts and tax table
Backhed_ARG_SSU_length_norm_res<-t(t(Backhed_arg_length_norm_res)/Backhed_SSU_counts[,1])*1541
Backhed_ARG_SSU_length_norm_res[is.na(Backhed_ARG_SSU_length_norm_res)]<-0
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)


#Read in sample data, sample data only from 98 samples
Backhed_sample_data<-read.csv("Backhed_newborn_metadata.csv", header = T, sep = ";")
row.names(Backhed_sample_data)<-Backhed_sample_data$Run

Backhed_sample_data$Antibiotics<-FALSE

Backhed_ARG<-phyloseq(otu_table(Backhed_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(Backhed_sample_data), tax_table(as.matrix(ARG_tax)))

Backhed_ARG_unnorm<-phyloseq(otu_table(Backhed_ARG_bt_res, taxa_are_rows = T), sample_data(Backhed_sample_data), tax_table(as.matrix(ARG_tax)))

Backhed_SSU_counts_unnorm<-Backhed_SSU_counts[(rownames(Backhed_SSU_counts)%in%c(sample_names(sample_data(Backhed_ARG_unnorm)))),]


#Using only resfinder hits remove genes which have no matches and make phyloseq object
Backhed_ARG_mod2_res<-otu_table(Backhed_ARG)[rowSums(otu_table(Backhed_ARG))>0]
Backhed_ARG_PHY_res<-phyloseq(Backhed_ARG_mod2_res, sample_data(Backhed_sample_data), tax_table(as.matrix(ARG_tax)))
#Backhed_ARG_PHY_res2<-(subset_taxa(Backhed_ARG_PHY_res, !V3%in%c("oqxB")))


#Add the explanatory variables to a dataframe

#Remove samples which have no info on feeding, 86 samples
Backhed_ARG_PHY_res2_diet<-subset_samples(Backhed_ARG_PHY_res, !feeding.practice.first.week%in%c(""))

Backhed_ARG_unnorm<-subset_samples(Backhed_ARG_unnorm, !feeding.practice.first.week%in%c(""))

#Remove samples which have been taken after the first week, 70
Backhed_ARG_PHY_res2_diet<-subset_samples(Backhed_ARG_PHY_res2_diet, Age.at.sample.Newborn..days.<8)
Backhed_ARG_unnorm<-subset_samples(Backhed_ARG_unnorm, Age.at.sample.Newborn..days.<8)

#Subset samples that have small libraries based on SSU counts and mapping

Backhed_ARG_PHY_res2_diet<-subset_samples(Backhed_ARG_PHY_res2_diet, Backhed_SSU_counts>1000)
Backhed_ARG_unnorm<-subset_samples(Backhed_ARG_unnorm, SSU_counts>1000)

Backhed_ARG_res2_diet<-otu_table(Backhed_ARG_PHY_res2_diet)[rowSums(otu_table(Backhed_ARG_PHY_res2_diet))>0]

Backhed_ARG_PHY_res2_diet<-phyloseq(otu_table(Backhed_ARG_res2_diet), sample_data(Backhed_sample_data), tax_table(as.matrix(ARG_tax)))

sample_data(Backhed_ARG_PHY_res2_diet)$Twin<-"n"


#Make dataframe

df_back <- data.frame(SUM=sample_sums(otu_table(Backhed_ARG_PHY_res2_diet)), Formula=sample_data(Backhed_ARG_PHY_res2_diet)$feeding.practice.first.week, Gest_age=40, Age=sample_data(Backhed_ARG_PHY_res2_diet)$Age.at.sample.Newborn..days., Delivery_mode=sample_data(Backhed_ARG_PHY_res2_diet)$Delivery.mode,
                      Country="Denmark",    Study="Backhed",         Antibiotics=sample_data(Backhed_ARG_PHY_res2_diet)$Antibiotics)


df_back$Formula<-!df_back$Formula%in%c("exclusively breastfeeding")

fit<-(glm(SUM~Age+Formula+Delivery_mode, data=df_back, family="Gamma"(link="log")))
summ(fit, exp = TRUE)



#Read in ARG mapping data for the late infancy samples and mothers

Backhed_ARG_bt_res<-as.matrix(read.table("Backhed_motinf_ARG_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F))


ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Backhed_ARG_bt_res<-Backhed_ARG_bt_res[row.names(Backhed_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Backhed_arg_length_norm_res <- Backhed_ARG_bt_res/ARG_lengths_res[,1]

Backhed_ARG_bt<-Backhed_ARG_bt_res[row.names(Backhed_ARG_bt_res)%in%row.names(ARG_lengths),]

#Read in SSU_counts
Backhed_SSU_counts_motinf <- as.matrix(read.table("Backhed_Mothers_Infants_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))

#Read in SSU counts and tax table
Backhed_ARG_SSU_length_norm_res<-t(t(Backhed_arg_length_norm_res)/Backhed_SSU_counts_motinf[,1])*1541
Backhed_ARG_SSU_length_norm_res[is.na(Backhed_ARG_SSU_length_norm_res)]<-0
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)

#Read in sample data from mothers and infants
Backhed_sample_data_motinf<-read.csv("Bäckhed_SraRunTable_motINf.txt", header = T, sep ="\t")
row.names(Backhed_sample_data)<-Backhed_sample_data$Run
Backhed_metadata<-read.csv("Backhed_metadata_BYSTUDYID.csv", header=T, sep=";")
Backhed_metadata <- Backhed_metadata %>% dplyr::rename(StudyID=STUDY.ID ) 

Backhed_motinf_sampledata<-merge(Backhed_sample_data_motinf, Backhed_metadata, by="StudyID")
row.names(Backhed_motinf_sampledata)<-Backhed_motinf_sampledata$Run

#Add ARG sum at birth to sample_data
sample_data(Backhed_ARG_PHY_res2_diet)$SUM=sample_sums(Backhed_ARG_PHY_res2_diet)

sample_data_temp<-data.frame(sample_data(Backhed_ARG_PHY_res2_diet))
temp<-base::merge(sample_data_temp, Backhed_motinf_sampledata, by.x="STUDY.ID", by.y="StudyID")

rownames(temp)<-temp$Run.y

Backhed_ARG<-phyloseq(otu_table(Backhed_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(temp), tax_table(as.matrix(ARG_tax)))

Backhed_ARG_4M<-subset_samples(Backhed_ARG, Run.y%in%c(row.names(sample_data(Backhed_ARG))[grep(x = sample_data(Backhed_ARG)$Alias2, pattern = "_4M")]))

Backhed_ARG_12M<-subset_samples(Backhed_ARG, Run.y%in%c(row.names(sample_data(Backhed_ARG))[grep(x = sample_data(Backhed_ARG)$Alias2, pattern = "_12M")]))

Backhed_ARG_M<-subset_samples(Backhed_ARG, Run.y%in%c(row.names(sample_data(Backhed_ARG))[grep(x = sample_data(Backhed_ARG)$Alias2, pattern = "_M")]))


df_4M<-data.frame(SUM=sample_data(Backhed_ARG_4M)$SUM, SUM12M=sample_sums(Backhed_ARG_12M), SUM4M=sample_sums(Backhed_ARG_4M), SUMM=sample_sums(Backhed_ARG_M), Formula=!sample_data(Backhed_ARG_12M)$feeding.practice.first.week.x=="exclusively breastfeeding", Formula4M=!sample_data(Backhed_ARG_12M)$feeding.practice.4M=="exclusively breastfeeding", BM12M=sample_data(Backhed_ARG_12M)$Any.breastfeeding.12.M=="any breastfeeding", Delivery_mode=sample_data(Backhed_ARG_12M)$Delivery.mode.y, AB_4M=sample_data(Backhed_ARG_12M)$Antibiotic.treatment.to.infant.0.4M...........times., AB_12M=sample_data(Backhed_ARG_12M)$Antibiotic.treatment.to.infant.4.12M.....times.,
                  RUN=sample_names(Backhed_ARG_12M), Age=sample_data(Backhed_ARG_12M)$Age.at.Sampling.12M.y)

fit_12M <-glm(SUM12M~SUMM+SUM+Age+SUM4M+AB_12M+Formula4M, data=df_4M, family = Gamma(link="log"))

summ(fit_12M, exp=TRUE, digits = 4)

```
# Banfield
https://doi.org/10.1128%2FmSystems.00123-17


```{r, results="hide", messages =FALSE, warnings =FALSE}
Banfield_ARG_bt_res<-as.matrix(read.table("ARG_genemat_banfield_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Banfield_ARG_bt_res<-Banfield_ARG_bt_res[row.names(Banfield_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Banfield_arg_length_norm_res <- Banfield_ARG_bt_res/ARG_lengths_res[,1]


#Read in SSU_counts, this is always >1000
Banfield_SSU_counts <- as.matrix(read.table("banfield_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))



#Read in SSU counts and tax table
Banfield_ARG_SSU_length_norm_res<-t(t(Banfield_arg_length_norm_res)/Banfield_SSU_counts[,1])*1541

Banfield_ARG_SSU_length_norm_res[is.na(Banfield_ARG_SSU_length_norm_res)]<-0

ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)


#Read in sample data
banfield_infant_metadata<-read.csv("/Users/kparnane/Documents/NEC/Analysis/banfield_infant_metadata_2.csv", sep = ";")

banfield_SRA_table<-read.table("/Users/kparnane/Documents/NEC/Analysis/SraRunTable_banfield.txt", sep = "\t", header=T)

banfield_SRA_table_cubs<-banfield_SRA_table[banfield_SRA_table$Run%in%colnames(Banfield_ARG_SSU_length_norm_res),]

banfield_infant_metadata_cubs<-banfield_infant_metadata[banfield_infant_metadata$infant%in%banfield_SRA_table_cubs$infant,]

sample_data_banfield<-merge(banfield_infant_metadata_cubs, data.frame(infant=banfield_SRA_table_cubs$infant, Run=banfield_SRA_table_cubs$Run, age=banfield_SRA_table_cubs$DOL), by="infant")

sample_data_banfield$Antibiotics<-!sample_data_banfield$antibiotics_taken=="[]"

row.names(sample_data_banfield)<-sample_data_banfield$Run



sample_data_banfield2<-read.table(file="Banfield_sampledata.txt", sep="\t", row.names = 1, header=TRUE)


sample_data_banfield2$Mat_AB<-sample_data_banfield2$maternal_ab

Banfield_ARG<-phyloseq(otu_table(Banfield_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(sample_data_banfield2), tax_table(as.matrix(ARG_tax)))



#Using only resfinder hits remove genes which have no matches and make phyloseq object
Banfield_ARG_mod2_res<-otu_table(Banfield_ARG)[rowSums(otu_table(Banfield_ARG))>0]
Banfield_ARG_PHY_res<-phyloseq(Banfield_ARG_mod2_res, sample_data(sample_data_banfield2), tax_table(as.matrix(ARG_tax)))
#Banfield_ARG_PHY_res2<-(subset_taxa(Banfield_ARG_PHY_res, !V3%in%c("oqxB")))
Banfield_ARG_PHY_res2<-Banfield_ARG_PHY_res
Banfield_ARG_unnorm<-phyloseq(otu_table(Banfield_ARG_bt_res, taxa_are_rows = T), sample_data(sample_data_banfield), tax_table(as.matrix(ARG_tax)))

sample_data(Banfield_ARG_unnorm)$SSU_counts<-Banfield_SSU_counts[,1]




df_ban <- data.frame(SUM=round(mean(Banfield_SSU_counts)*(sample_sums(otu_table(Banfield_ARG_PHY_res2)))),
                 Formula=sample_data(Banfield_ARG_PHY_res2)$feeding%in%c("Formula", "Combination"),
                 Gest_age=sample_data(Banfield_ARG_PHY_res2)$gestational_age,
                 Age=sample_data(Banfield_ARG_PHY_res2)$age,
                 Delivery_mode=sample_data(Banfield_ARG_PHY_res2)$birth_mode,
                 Country="US",
                 mat_ab=sample_data(Banfield_ARG_PHY_res2)$maternal_ab,
                 Study="Banfield",
                 NEC=sample_data(Banfield_ARG_PHY_res2)$NEC,
                 Antibiotics=sample_data(Banfield_ARG_PHY_res2)$Antibiotics,
                 antibiotics_taken=sample_data(Banfield_ARG_PHY_res2)$antibiotics_taken,
                 sepsis=sample_data(Banfield_ARG_PHY_res2)$sepsis,
                 Pair=sample_data(Banfield_ARG_PHY_res2)$Pair,
                 diseases=!sample_data(Banfield_ARG_PHY_res2)$diseases=="[]")

fit<-(glm(SUM~mat_ab+Delivery_mode+sepsis+diseases+Antibiotics+Age+Gest_age+Formula, data=df_ban, family="Gamma"(link="log")))

fit<-(glm(SUM~Age+Gest_age+Formula, data=df_ban, family="Gamma"(link="log")))
summ(fit, exp = TRUE)




```
# Dsouza, Dantas lab Nature medicine
(Second author last name used for shortness for this publication in the code block. First author is Baumann-Dudenhoeffer, A.M.)
https://doi.org/10.1038/s41591-018-0216-2
```{r, results="hide", messages =FALSE, warnings =FALSE}

# Read in ARG data
Dsouza_ARG_bt_res<-as.matrix(read.table("Dsouza_ARG_resfinder_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

# Read in 16S counts

#Read in SSU_counts
Dsouza_SSU_counts <- as.matrix(read.table("Dsouza_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Dsouza_ARG_bt_res<-Dsouza_ARG_bt_res[row.names(Dsouza_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Dsouza_arg_length_norm_res <- Dsouza_ARG_bt_res/ARG_lengths_res[,1]

# Read in sample data
Dsouza_sample_data<-read.table("/Users/kparnane/Documents/NEC/Analysis/SraRunTable_Dsouza.txt", row.names = 1, header=T, sep="\t")
Dsouza_sample_data$Antibiotics<-FALSE
Dsouza_sample_data$Mat_AB<-!Dsouza_sample_data$host_maternal_peripartum_abx==""


#Read in SSU counts and tax table
Dsouza_ARG_SSU_length_norm_res<-t(t(Dsouza_arg_length_norm_res)/Dsouza_SSU_counts[,1])*1541

Dsouza_ARG_SSU_length_norm_res[is.na(Dsouza_ARG_SSU_length_norm_res)]<-0

ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)


#Make phyloseq
Dsouza_ARG_unnorm<-phyloseq(otu_table(Dsouza_ARG_bt_res, taxa_are_rows = T), sample_data(Dsouza_sample_data), tax_table(as.matrix(ARG_tax)))

Dsouza_ARG_res<-phyloseq(otu_table(Dsouza_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(Dsouza_sample_data), tax_table(as.matrix(ARG_tax)))
Dsouza_ARG_PHY_res2<-Dsouza_ARG_res

#Remove samples with small libraries

Dsouza_ARG_PHY_res2<-subset_samples(Dsouza_ARG_PHY_res2, (Dsouza_SSU_counts>1000))
sample_data(Dsouza_ARG_unnorm)$SSU_counts<-Dsouza_SSU_counts[,1]
Dsouza_ARG_unnorm<-subset_samples(Dsouza_ARG_unnorm, (Dsouza_SSU_counts>1000))

Dsouza_ARG_PHY_res2<-subset_samples(Dsouza_ARG_PHY_res2, sample_sums(Dsouza_ARG_PHY_res2)>0)

sample_data(Dsouza_ARG_PHY_res2)$Pair<-gsub(sample_data(Dsouza_ARG_PHY_res2)$Sample_Name, pattern = "._0", replacement = "")

#Make df for GLM

df_dsouza<-data.frame(SUM=(sample_sums(Dsouza_ARG_unnorm)), 
                 Gest_age=sample_data(Dsouza_ARG_PHY_res2)$host_gestationalage_atdelivery,
                  Formula=sample_data(Dsouza_ARG_PHY_res2)$Formula,
                 Age=sample_data(Dsouza_ARG_PHY_res2)$host_age_days_at_time_of_survey,
                 Country="US",
                 Gender=sample_data(Dsouza_ARG_PHY_res2)$host_sex,
                 Study="Dsouza",
                 Delivery_mode=sample_data(Dsouza_ARG_PHY_res2)$host_del_route,
                 Antibiotics=sample_data(Dsouza_ARG_PHY_res2)$Antibiotics,
                 Mat_ab=sample_data(Dsouza_ARG_PHY_res2)$Mat_AB)

fit<-(glm(SUM~Gest_age+Age+Mat_ab+Formula, data=df_dsouza, family="Gamma"(link="log")))
summ(fit, exp = TRUE)


#Infancy

# Read in ARG data
Dsouza_ARG_bt_res<-as.matrix(read.table("Dsouza_inf_ARG_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

# Read in 16S counts

#Read in SSU_counts
Dsouza_SSU_counts_inf <- as.matrix(read.table("Dsouza_Infants_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Dsouza_ARG_bt_res<-Dsouza_ARG_bt_res[row.names(Dsouza_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Dsouza_arg_length_norm_res <- Dsouza_ARG_bt_res/ARG_lengths_res[,1]

# Read in sample data
Dsouza_sample_data<-read.table("/Users/kparnane/Documents/NEC/Analysis/Dsouza_infancy_sampledata.txt", row.names = 1, header=T, sep="\t")
Dsouza_sample_data$Mat_AB<-!Dsouza_sample_data$host_maternal_peripartum_abx==""


#Read in SSU counts and tax table
Dsouza_ARG_SSU_length_norm_res<-t(t(Dsouza_arg_length_norm_res)/Dsouza_SSU_counts_inf[,1])*1541

Dsouza_ARG_SSU_length_norm_res[is.na(Dsouza_ARG_SSU_length_norm_res)]<-0

ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)


#Make phyloseq

Dsouza_ARG<-phyloseq(otu_table(Dsouza_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(Dsouza_sample_data), tax_table(as.matrix(ARG_tax)))

#Merge with newborn data
Dsouza_all<-merge_phyloseq(Dsouza_ARG_PHY_res2, Dsouza_ARG)

#Subset samples
sample_data(Dsouza_all)$SUM<-sample_sums(Dsouza_all)
Dsouza_NB<-subset_samples(Dsouza_all, host_age_days_at_time_of_survey<36)
Dsouza_4M<-subset_samples(Dsouza_ARG, host_age_days_at_time_of_survey>100&host_age_days_at_time_of_survey<200)
sample_data(Dsouza_4M)$SUM<-sample_sums(Dsouza_4M)
sample_data(Dsouza_4M)$Host_diet4m<-sample_data(Dsouza_4M)$Host_Diet
Dsouza_8M<-subset_samples(Dsouza_ARG, host_age_days_at_time_of_survey>200)
sample_data(Dsouza_8M)$SUM<-sample_sums(Dsouza_8M)
sample_data(Dsouza_8M)$Host_diet8m<-sample_data(Dsouza_8M)$Host_Diet

temp<-base::merge(data.frame(sample_data(Dsouza_NB)), data.frame(sample_data(Dsouza_4M)), by="host_subject_id")

temp2<-temp %>% rename(SUMNB=SUM.x, SUM4M=SUM.y)

temp3<-base::merge(temp2, data.frame(sample_data(Dsouza_8M)), by="host_subject_id")

df_Dsouza_inf<-data.frame(SUM=temp3$SUMNB, SUM4M=temp3$SUM4M, Age=temp3$host_age_days_at_time_of_survey.x, gest_Age=temp3$host_gestationalage_atdelivery.x, SUM12M=temp3$SUM, FormulaNB=temp3$Formula, Formula4M=!temp3$Host_diet4m=="Exclusively Breastfed\\, N/A", Formula8M=!temp3$Host_diet8m=="Exclusively Breastfed\\, N/A")

df_joined<-(full_join(df_Dsouza_inf, df_4M))

summ(glm(SUM4M~Age+Formula4M+SUM, data=df_joined, family = Gamma(link="log")), digits=4, exp=TRUE)

summ(glm(SUM12M~Formula4M+SUM4M, data=df_joined, family = Gamma(link="log")), digits=4, exp=TRUE)

```
# Gibson

https://doi.org/10.1038/nmicrobiol.2016.24

```{r, results="hide", messages =FALSE, warnings =FALSE}
# Read in ARG data
Gibson_ARG_bt_res<-as.matrix(read.table("Gibson_ARG_resfinder_genemat.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

# Read in 16S counts

#Read in SSU_counts
Gibson_SSU_counts <- as.matrix(read.table("Gibson_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Gibson_ARG_bt_res<-Gibson_ARG_bt_res[row.names(Gibson_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Gibson_arg_length_norm_res <- Gibson_ARG_bt_res/ARG_lengths_res[,1]

# Read in sample data
Gibson_sample_data<-read.csv("/Users/kparnane/Documents/NEC/Analysis/Gibson_metadata_foranalysis.csv", row.names = 1, header=T, sep=";")

Gibson_sample_data$Mat_AB<-rowSums(Gibson_sample_data[,47:57])>0

Gibson_sample_data$Antibiotics<-Gibson_sample_data$Total_abx>0


#Read in SSU counts and tax table
Gibson_ARG_SSU_length_norm_res<-t(t(Gibson_arg_length_norm_res)/Gibson_SSU_counts[,1])*1541

Gibson_ARG_SSU_length_norm_res[is.na(Gibson_ARG_SSU_length_norm_res)]<-0

ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
#Modify tax table
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)


#Make phyloseq
Gibson_ARG_res<-phyloseq(otu_table(Gibson_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(Gibson_sample_data), tax_table(as.matrix(ARG_tax)))
Gibson_ARG_PHY_res2<-Gibson_ARG_res
#Gibson_ARG_PHY_res2<-(subset_taxa(Gibson_ARG_res, !V3%in%c("oqxB")))

sample_data(Gibson_ARG_PHY_res2)$Mat_AB<-sample_data(Gibson_ARG_PHY_res2)$Mat_AB>0
Gibson_ARG_unnorm<-phyloseq(otu_table(Gibson_ARG_bt_res, taxa_are_rows = T), sample_data(Gibson_sample_data), tax_table(as.matrix(ARG_tax)))

sample_data(Gibson_ARG_PHY_res2)$SSU_counts<-Gibson_SSU_counts[,1]



#Remove samples with small libraries

Gibson_ARG_PHY_res2<-subset_samples(Gibson_ARG_PHY_res2, (Gibson_SSU_counts>1000))

Gibson_ARG_unnorm<-subset_samples(Gibson_ARG_unnorm, (Gibson_SSU_counts>1000))

Gibson_ARG_PHY_res2<-subset_samples(Gibson_ARG_PHY_res2, sample_sums(Gibson_ARG_PHY_res2)>0)

sample_data(Gibson_ARG_PHY_res2)$Pair<-gsub(sample_data(Gibson_ARG_PHY_res2)$Individual, pattern = "\\.0.", replacement="")



#Make df for GLM

subs_Gibson<-subset_samples(Gibson_ARG_PHY_res2)

df_gibson<-data.frame(SUM=sample_sums(subs_Gibson), 
                 Gest_age=sample_data(subs_Gibson)$Gestational_Age,
                Formula=sample_data(subs_Gibson)$Formula>0,
                 Age=sample_data(subs_Gibson)$Day_of_Life,
                Delivery_mode=sample_data(subs_Gibson)$Delivery_Mode,
                 Country="US",
                Study="Gibson",
                Gender=sample_data(subs_Gibson)$Gender,
                Antibiotics=sample_data(subs_Gibson)$Antibiotics,
                Individual=sample_data(subs_Gibson)$Individual,
                NEC=sample_data(subs_Gibson)$necrotizing_enterocolitis_status,
                SSU_counts=sample_data(subs_Gibson)$SSU_counts,
                Mat_AB=sample_data(subs_Gibson)$Mat_AB)

#write.csv(df_gibson, "gibson_data.csv", quote=FALSE)
#df_gibson<-read.csv("gibson_data.csv", sep=";")


fit<-(glm(SUM~SSU_counts+Mat_AB+Age+Antibiotics+Gest_age+Formula, data=df_gibson, family="Gamma"(link="log")))
summ(fit, exp=TRUE)


```

# Gather data for phyloseq for meta-analysis

The following code is for harmonizing and compiling metadata for meta-analysis steps.

```{r, results="hide", messages =FALSE, warnings =FALSE}
#Add SSU counts to metadata table
ARG_PHY_res2<-ARG_PHY_res
sample_data(ARG_PHY_res2)$SSU_counts<-SSU_counts[,1]
Backhed_SSU_counts2<-Backhed_SSU_counts[(rownames(Backhed_SSU_counts)%in%c(sample_names(sample_data(Backhed_ARG_PHY_res2_diet)))),]
sample_data(Backhed_ARG_PHY_res2_diet)$SSU_counts<-Backhed_SSU_counts2

Gibson_SSU_counts2<-Gibson_SSU_counts[(rownames(Gibson_SSU_counts)%in%c(sample_names(sample_data(Gibson_ARG_PHY_res2)))),]
sample_data(Gibson_ARG_PHY_res2)$SSU_counts<-Gibson_SSU_counts2

sample_data(Banfield_ARG_PHY_res2)$SSU_counts<-Banfield_SSU_counts[,1]

Dsouza_SSU_counts2<-Dsouza_SSU_counts[(rownames(Dsouza_SSU_counts)%in%c(sample_names(sample_data(Dsouza_ARG_PHY_res2)))),]
sample_data(Dsouza_ARG_PHY_res2)$SSU_counts<-Dsouza_SSU_counts2

#Add study to metadata tables
sample_data(ARG_PHY_res2)$Study="NEC"
sample_data(Backhed_ARG_PHY_res2_diet)$Study="Backhed"
sample_data(Gibson_ARG_PHY_res2)$Study="Gibson"
sample_data(Banfield_ARG_PHY_res2)$Study="Rahman"
sample_data(Dsouza_ARG_PHY_res2)$Study="Baumann-Dudenhoeffer"


#Add formula to metadatatables
sample_data(Backhed_ARG_PHY_res2_diet)$Formula2<-!sample_data(Backhed_ARG_PHY_res2_diet)$feeding.practice.first.week=="exclusively breastfeeding"

sample_data(Gibson_ARG_PHY_res2)$Formula2<-sample_data(Gibson_ARG_PHY_res2)$Formula>0

sample_data(Dsouza_ARG_PHY_res2)$Formula2<-sample_data(Dsouza_ARG_PHY_res2)$Formula

sample_data(ARG_PHY_res2)$Formula2<-sample_data(ARG_PHY_res2)$Formula=="y"
sample_data(Banfield_ARG_PHY_res2)$Formula2<-!sample_data(Banfield_ARG_PHY_res2)$feeding=="Breast"

#Add maternal AB to metadata tables
sample_data(Backhed_ARG_PHY_res2_diet)$Mat_AB<-"NA"
sample_data(ARG_PHY_res2)$Mat_AB<-sample_data(ARG_PHY_res2)$Mat_AB=="Y"
sample_data(Dsouza_ARG_PHY_res2)$Mat_AB
sample_data(Gibson_ARG_PHY_res2)$Mat_AB
sample_data(Banfield_ARG_PHY_res2)$Mat_AB

#Add gestational age to metadata
sample_data(Dsouza_ARG_PHY_res2)$Gest_age <-
sample_data(Dsouza_ARG_PHY_res2)$host_gestationalage_atdelivery

sample_data(Banfield_ARG_PHY_res2)$Gest_age <-
sample_data(Banfield_ARG_PHY_res2)$gestational_age

sample_data(Gibson_ARG_PHY_res2)$Gest_age <- sample_data(Gibson_ARG_PHY_res2)$Gestational_Age

sample_data(Backhed_ARG_PHY_res2_diet)$Gest_age<-40

#Add age and twin to metadata
sample_data(Banfield_ARG_PHY_res2)$Age <- sample_data(Banfield_ARG_PHY_res2)$age

sample_data(Banfield_ARG_PHY_res2)$Twin<-as.character(sample_data(Banfield_ARG_PHY_res2)$Twin)

sample_data(Dsouza_ARG_PHY_res2)$Twin<-as.character(sample_data(Dsouza_ARG_PHY_res2)$Twin)


sample_data(Gibson_ARG_PHY_res2)$Twin<-as.character(sample_data(Gibson_ARG_PHY_res2)$Twin)

sample_data(Dsouza_ARG_PHY_res2)$Age <- sample_data(Dsouza_ARG_PHY_res2)$host_age_days_at_time_of_survey

sample_data(Backhed_ARG_PHY_res2_diet)$Age<-sample_data(Backhed_ARG_PHY_res2_diet)$Age.at.sample.Newborn..days.

sample_data(ARG_PHY_res2)$Age<-sample_data(ARG_PHY_res2)$Age_sampling
  
sample_data(Gibson_ARG_PHY_res2)$Age<-sample_data(Gibson_ARG_PHY_res2)$Day_of_Life

#Add delivery method
sample_data(ARG_PHY_res2)$Delivery_mode<-as.character(sample_data(ARG_PHY_res2)$Delivery.Method)
sample_data(Banfield_ARG_PHY_res2)$Delivery_mode <- as.character(sample_data(Banfield_ARG_PHY_res2)$birth_mode)
sample_data(Backhed_ARG_PHY_res2_diet)$Delivery_mode<-as.character(sample_data(Backhed_ARG_PHY_res2_diet)$Delivery.mode)
sample_data(Dsouza_ARG_PHY_res2)$Delivery_mode<-as.character(sample_data(Dsouza_ARG_PHY_res2)$host_del_route)
sample_data(Gibson_ARG_PHY_res2)$Delivery_mode<-as.character(sample_data(Gibson_ARG_PHY_res2)$Delivery_Mode)

#Add family
sample_data(ARG_PHY_res2)$Pair
sample_data(Banfield_ARG_PHY_res2)


#Add host gender when missing
sample_data(Banfield_ARG_PHY_res2)$Gender<-sample_data(Banfield_ARG_PHY_res2)$gender
sample_data(Backhed_ARG_PHY_res2_diet)$Gender<-sample_data(Backhed_ARG_PHY_res2_diet)$GENDER
sample_data(Backhed_ARG_PHY_res2_diet)$Gender<-gsub("boy", "M", sample_data(Backhed_ARG_PHY_res2_diet)$Gender)
sample_data(Backhed_ARG_PHY_res2_diet)$Gender<-gsub("girl", "F", sample_data(Backhed_ARG_PHY_res2_diet)$Gender)
sample_data(Dsouza_ARG_PHY_res2)$Gender<-sample_data(Dsouza_ARG_PHY_res2)$host_sex
sample_data(Dsouza_ARG_PHY_res2)$Gender<-gsub("female", "F", sample_data(Dsouza_ARG_PHY_res2)$Gender)
sample_data(Dsouza_ARG_PHY_res2)$Gender<-gsub("male", "M", sample_data(Dsouza_ARG_PHY_res2)$Gender)
sample_data(Gibson_ARG_PHY_res2)$Gender<-as.character(sample_data(Gibson_ARG_PHY_res2)$Gender)
sample_data(Banfield_ARG_PHY_res2)$Gender<-as.character(sample_data(Banfield_ARG_PHY_res2)$Gender)

#Add maternal pre-eclampsia
sample_data(ARG_PHY_res2)$Maternal_preeclamp<-gsub("y", TRUE, sample_data(ARG_PHY_res2)$Maternal_preeclamp)
sample_data(ARG_PHY_res2)$Maternal_preeclamp<-gsub("n", FALSE, sample_data(ARG_PHY_res2)$Maternal_preeclamp)
sample_data(Banfield_ARG_PHY_res2)$Maternal_preeclamp<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$Maternal_preeclamp<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$Maternal_preeclamp<-sample_data(Dsouza_ARG_PHY_res2)$host_maternal_preeclampsia
sample_data(Gibson_ARG_PHY_res2)$Maternal_preeclamp<-sample_data(Gibson_ARG_PHY_res2)$pre.eclampsia
sample_data(Gibson_ARG_PHY_res2)$Maternal_preeclamp<-gsub("Yes", TRUE, sample_data(Gibson_ARG_PHY_res2)$Maternal_preeclamp)
sample_data(Gibson_ARG_PHY_res2)$Maternal_preeclamp<-gsub("No", FALSE, sample_data(Gibson_ARG_PHY_res2)$Maternal_preeclamp)

#Add fortifier
sample_data(ARG_PHY_res2)$Fortifier<-sample_data(ARG_PHY_res2)$Milk_Fortifier=="y"
sample_data(Gibson_ARG_PHY_res2)$Fortifier<-sample_data(Gibson_ARG_PHY_res2)$Fortification>0




#Add SUM

sample_data(ARG_PHY_res2)$SUM<-round(sample_sums(ARG_PHY_res2))
sample_data(Banfield_ARG_PHY_res2)$SUM<-round(sample_sums(Banfield_ARG_PHY_res2))
sample_data(Dsouza_ARG_PHY_res2)$SUM<-round(sample_sums(Dsouza_ARG_PHY_res2))
sample_data(Gibson_ARG_PHY_res2)$SUM<-round(sample_sums(Gibson_ARG_PHY_res2))
sample_data(Backhed_ARG_PHY_res2_diet)$SUM<-round(sample_sums(Backhed_ARG_PHY_res2_diet))

#Add country and sample when missing
sample_data(ARG_PHY_res2)$Country<-"US"
sample_data(Banfield_ARG_PHY_res2)$Country<-"US"
sample_data(Banfield_ARG_PHY_res2)$Sample<-sample_data(Banfield_ARG_PHY_res2)$Run
sample_data(Dsouza_ARG_PHY_res2)$Country<-"Luxemburg"
sample_data(Dsouza_ARG_PHY_res2)$Sample<-row.names(sample_data(Dsouza_ARG_PHY_res2))
sample_data(Backhed_ARG_PHY_res2_diet)$Country<-"Denmark"
sample_data(Backhed_ARG_PHY_res2_diet)$Sample<-sample_data(Backhed_ARG_PHY_res2_diet)$Run
sample_data(Gibson_ARG_PHY_res2)$Country<-"US"
sample_data(Gibson_ARG_PHY_res2)$Sample<-row.names(sample_data(Gibson_ARG_PHY_res2))

#Add maternal antibiotics
sample_data(ARG_PHY_res2)$Mat_AB<-sample_data(ARG_PHY_res2)$Mat_AB=="Y"
sample_data(Banfield_ARG_PHY_res2)$Mat_AB
sample_data(Gibson_ARG_PHY_res2)$Mat_AB
sample_data(Dsouza_ARG_PHY_res2)$Mat_AB
sample_data(Backhed_ARG_PHY_res2_diet)$Mat_AB<-NA

#Add NEC
sample_data(ARG_PHY_res2)$NEC<-FALSE
sample_data(Gibson_ARG_PHY_res2)$NEC<-sample_data(Gibson_ARG_PHY_res2)$necrotizing_enterocolitis_status=="Yes"
sample_data(Backhed_ARG_PHY_res2_diet)$NEC<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$NEC<-FALSE
sample_data(Banfield_ARG_PHY_res2)[grep("NEC", sample_data(Banfield_ARG_PHY_res2)$diseases),]$NEC<-TRUE
sample_data(Banfield_ARG_PHY_res2)[grep("NEC", sample_data(Banfield_ARG_PHY_res2)$diseases, invert = TRUE),]$NEC<-FALSE


#Add antibiotic types when missing

sample_data(ARG_PHY_res2)$cefotaxime<-FALSE
sample_data(ARG_PHY_res2)$clindamycin<-FALSE
sample_data(ARG_PHY_res2)$ofloxacin<-FALSE
sample_data(ARG_PHY_res2)$unknown<-FALSE
sample_data(ARG_PHY_res2)$cefazolin<-FALSE
sample_data(ARG_PHY_res2)$meropenem<-FALSE
sample_data(ARG_PHY_res2)$trimethoprim_sulfamathoxazole<-FALSE
sample_data(ARG_PHY_res2)$ticarcillin_clavulanate<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$gentamycin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$ampicillin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$vancomycin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$cefotaxime<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$clindamycin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$ofloxacin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$meropenem<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$unknown<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$cefazolin<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$ticarcillin_clavulanate<-FALSE
sample_data(Dsouza_ARG_PHY_res2)$trimethoprim_sulfamathoxazole<-FALSE

sample_data(Backhed_ARG_PHY_res2_diet)$trimethoprim_sulfamathoxazole<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$gentamycin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$ampicillin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$vancomycin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$cefotaxime<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$clindamycin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$ofloxacin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$meropenem<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$unknown<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$cefazolin<-FALSE
sample_data(Backhed_ARG_PHY_res2_diet)$ticarcillin_clavulanate<-FALSE

sample_data(Gibson_ARG_PHY_res2)$gentamycin<-sample_data(Gibson_ARG_PHY_res2)$Gentamicin>0
sample_data(Gibson_ARG_PHY_res2)$ampicillin<-sample_data(Gibson_ARG_PHY_res2)$Ampicillin>0
sample_data(Gibson_ARG_PHY_res2)$vancomycin<-sample_data(Gibson_ARG_PHY_res2)$Vancomycin>0
sample_data(Gibson_ARG_PHY_res2)$ticarcillin_clavulanate<-sample_data(Gibson_ARG_PHY_res2)$Ticarcillin.Clavulanate>0
sample_data(Gibson_ARG_PHY_res2)$clindamycin<-sample_data(Gibson_ARG_PHY_res2)$Clindamycin>0
sample_data(Gibson_ARG_PHY_res2)$cefotaxime<-sample_data(Gibson_ARG_PHY_res2)$Cefotaxime>0
sample_data(Gibson_ARG_PHY_res2)$cefazolin<-sample_data(Gibson_ARG_PHY_res2)$Cefazolin>0
sample_data(Gibson_ARG_PHY_res2)$ofloxacin<-FALSE
sample_data(Gibson_ARG_PHY_res2)$unknown<-FALSE
sample_data(Gibson_ARG_PHY_res2)$trimethoprim_sulfamathoxazole<-sample_data(Gibson_ARG_PHY_res2)$Trimethoprim.Sulfamathoxazole>0
sample_data(Gibson_ARG_PHY_res2)$meropenem<-sample_data(Gibson_ARG_PHY_res2)$Meropenem>0
sample_data(Banfield_ARG_PHY_res2)$meropenem<-FALSE
sample_data(Banfield_ARG_PHY_res2)$trimethoprim_sulfamathoxazole<-FALSE
sample_data(Banfield_ARG_PHY_res2)$ticarcillin_clavulanate<-FALSE
sample_data(Banfield_ARG_PHY_res2)$cefazolin<-FALSE
sample_data(Gibson_ARG_PHY_res2)$unkown<-FALSE

#Merge phyloseq
merged_phyloseq<-merge_phyloseq(ARG_PHY_res2, Banfield_ARG_PHY_res2, Dsouza_ARG_PHY_res2, Gibson_ARG_PHY_res2, Backhed_ARG_PHY_res2_diet)



sample_data(merged_phyloseq)$Fortifier<-as.character(sample_data(merged_phyloseq)$Fortifier)
sample_data(merged_phyloseq)$Fortifier<-ifelse(is.na(sample_data(merged_phyloseq)$Fortifier), 
             'FALSE', sample_data(merged_phyloseq)$Fortifier)


sample_data(merged_phyloseq)$Delivery_mode<-gsub("c-section", "CS", sample_data(merged_phyloseq)$Delivery_mode)
sample_data(merged_phyloseq)$Delivery_mode<-gsub("C-section", "CS", sample_data(merged_phyloseq)$Delivery_mode)
sample_data(merged_phyloseq)$Delivery_mode<-gsub("section", "CS", sample_data(merged_phyloseq)$Delivery_mode)
sample_data(merged_phyloseq)$Delivery_mode<-gsub("Cesarean", "CS", sample_data(merged_phyloseq)$Delivery_mode)
sample_data(merged_phyloseq)$Delivery_mode<-gsub("cesarean", "CS", sample_data(merged_phyloseq)$Delivery_mode) 
sample_data(merged_phyloseq)$Delivery_mode<-gsub("vaginal", "V", sample_data(merged_phyloseq)$Delivery_mode)
sample_data(merged_phyloseq)$Delivery_mode<-gsub("Vaginal", "V", sample_data(merged_phyloseq)$Delivery_mode)


#Keep only useful metadata
sample_data(merged_phyloseq)<-sample_data(merged_phyloseq)[,colnames(sample_data(merged_phyloseq))%in%c("Formula2", "Fortifier", "Delivery_mode", "Age", "Antibiotics", "Gest_age", "SSU_counts","gentamycin", "ampicillin", "vancomycin", "meropenem", "ticarcillin_clavulanate","ofloxacin", "NEC", "unknown", "cefazolin", "cefotaxime", "clindamycin", "Study", "Donor_Milk", "Sample", "Gender", "Country", "Twin", "Mat_AB", "Maternal_preeclamp", "Pair")]


otu_table(merged_phyloseq)<-otu_table(merged_phyloseq)[rowSums(is.na(otu_table(merged_phyloseq))) == 0,]

merged_phyloseq<-subset_taxa(merged_phyloseq, rowSums(otu_table(merged_phyloseq))>0)

merged_phyloseq<-subset_taxa(merged_phyloseq, !is.na(rowSums(otu_table(merged_phyloseq))))

write.table(sample_data(merged_phyloseq), file="Supplementary tables/Supplementary table 4-Meta-analysis sample data.txt", quote=FALSE, row.names = TRUE)


df<-data.frame(sample_data(merged_phyloseq))

df_F<-df[df$Formula2=="TRUE",]
df_NF<-df[df$Formula2=="FALSE",]


table(df_F$Delivery_mode)
summary(df_NF$Gest_age)

# Save RDS file
#saveRDS(merged_phyloseq, file="merged_phyloseq.rds")
```

# Meta-analysis analyses

```{r, results="hide", messages =FALSE, warnings =FALSE}

# Read RDS file for running the following code
readRDS(merged_phyloseq, file="merged_phyloseq.rds")

#Make data frame

df_meta<-data.frame(SUM=sample_sums(merged_phyloseq),
                    Gender=sample_data(merged_phyloseq)$Gender,
                    NEC=sample_data(merged_phyloseq)$NEC,
                    Formula=as.factor(sample_data(merged_phyloseq)$Formula2), Antibiotics=sample_data(merged_phyloseq)$Antibiotics, Delivery_mode=sample_data(merged_phyloseq)$Delivery_mode,
                    Fortifier=sample_data(merged_phyloseq)$Fortifier,
                    Age=sample_data(merged_phyloseq)$Age,
                    Gest_age=sample_data(merged_phyloseq)$Gest_age,
                    Study=sample_data(merged_phyloseq)$Study,
                    Sample=sample_names(merged_phyloseq),
                    Donor=sample_data(merged_phyloseq)$Donor_Milk,
                    Country=sample_data(merged_phyloseq)$Country,
                    SSU_counts=sample_data(merged_phyloseq)$SSU_counts,
                    Twin=sample_data(merged_phyloseq)$Twin,
                    Mat_AB=sample_data(merged_phyloseq)$Mat_AB,
                    gent=sample_data(merged_phyloseq)$gentamycin,
                    amp=sample_data(merged_phyloseq)$ampicillin,
                    van=sample_data(merged_phyloseq)$vancomycin,
                    meropenem=sample_data(merged_phyloseq)$meropenem,
                    cefotaxime=sample_data(merged_phyloseq)$cefotaxime,
                    ofloxacin=sample_data(merged_phyloseq)$ofloxacin,
                    clindamycin=sample_data(merged_phyloseq)$clindamycin,
                    cefazolin=sample_data(merged_phyloseq)$cefazolin,
                    Maternal_preeclamp=sample_data(merged_phyloseq)$Maternal_preeclamp)

#Save df for fortifier analysis
df_meta_fort<-df_meta[df_meta$Study%in%c("NEC", "Gibson"),]

df_meta_fort<-df_meta_fort[df_meta_fort$Formula=="FALSE",]

#Exclude fortifier samples that don't have formula as well
df_meta<-df_meta[!(df_meta$Fortifier=="TRUE"&df_meta$Formula=="FALSE"),]
```


# GLMs for meta-analysis

Generalized linear models for the meta-analysis dataset

```{r, results="hide", messages =FALSE, warnings =FALSE}

#Build models
df_meta$Formula<-as.logical(df_meta$Formula)
df_meta$Fortifier<-as.logical(df_meta$Fortifier)

fit_meta1<-glm(SUM~Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta1, exp=TRUE, digits=4)
summary(fit_meta1)

fit_meta2<-glm(SUM~Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta2, exp=TRUE, digits=4)
summary(fit_meta2)

fit_meta3<-glm(SUM~SSU_counts+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta3, exp=TRUE, digits=4)
summary(fit_meta3)

fit_meta4<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta4, exp=TRUE, digits = 4)
summary(fit_meta4)

fit_meta4_van<-glm(SUM~SSU_counts+Age+van+Gest_age+Fortifier+Formula, data=df_meta, family="Gamma"(link="log"))

summary(fit_meta4_van, exp=TRUE, digits = 4)
summ(fit_meta4_van, exp=TRUE ,digits=4)

fit_meta4_van_gender<-glm(SUM~SSU_counts+Age+van+Gender+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))


anova(fit_meta4, fit_meta4_van, fit_meta4_van_gender, test="Chisq")

library(lme4)
summary(glmer(SUM~Age+Gest_age+Formula + (1|Study), data=df_meta, family="Gamma"(link="log")))


fit_fullModel<-glm(SUM~Study+SSU_counts+Antibiotics+Mat_AB+Maternal_preeclamp+NEC+Delivery_mode+amp+van+meropenem+cefotaxime+clindamycin+ofloxacin+Age+Gest_age+Fortifier+Formula, data=df_meta, family="Gamma"(link="log"))

summary(fit_fullModel, exp=TRUE, digits=4)

(cor(model.matrix(fit_fullModel)[,-1]))

fit_meta4_AB<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta[df_meta$Antibiotics==TRUE,], family="Gamma"(link="log"))
summ(fit_meta4_AB, exp=TRUE)
summary(fit_meta4_AB)

fit_meta5<-glm(SUM~SSU_counts+Age+Gest_age+Antibiotics+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta5, exp=TRUE)
summary(fit_meta5)

anova(fit_meta4, fit_meta5, test="Chisq")

df_meta_pret_noab<-(df_meta[df_meta$Antibiotics=="FALSE"&df_meta$Gest_age<37,])
fit_meta4_pret_noab<-glm(SUM~SSU_counts+Age+Gest_age+Mat_AB+Antibiotics+Formula, data=df_meta_pret_noab, family="Gamma"(link="log"))
summary(fit_meta4_pret_noab)

fit_meta4_preterm<-glm(SUM~Formula+SSU_counts+Age+Gest_age, data=df_meta[df_meta$Gest_age<37,], family="Gamma"(link="log"))
summ(fit_meta4_preterm, exp=TRUE)
summary(fit_meta4_preterm)

fit_meta4_fullterm<-glm(SUM~Formula+SSU_counts+Delivery_mode+Age+Gest_age, data=df_meta[df_meta$Gest_age>=37,], family="Gamma"(link="log"))
summ(fit_meta4_fullterm, exp=TRUE)
summary(fit_meta4_fullterm)

fit_meta4_fullterm<-glm(SUM~Formula+Age, data=df_meta[df_meta$Gest_age>=37,], family="Gamma"(link="log"))
summ(fit_meta4_fullterm, exp=TRUE)
summary(fit_meta4_fullterm)


modelEffectSizes(fit_meta4_preterm)
modelPower(pc=3, pa=4, peta2 = 0.03,power=0.5)


pchisq(summary(fit_meta4)$deviance, 
           summary(fit_meta4)$df.residual
           )
#######################
#Plot predicted versus actual values
temp<-cbind(df_meta, 
      Mean = predict(fit_meta4, newdata=df_meta, type="response"), 
      SE = predict(fit_meta4, newdata=df_meta, type="response", se.fit=T)$se.fit
      )

fam<-family(fit_meta4)
str(fam)

ilink <- fam$linkinv
ilink

df_meta <- bind_cols(df_meta, setNames(as_tibble(predict(fit_meta4, df_meta, se.fit = TRUE)[1:2]), c('fit_link','se_link')))

df_meta <- mutate(df_meta, fit_resp = ilink(fit_link), right_upr = ilink(fit_link + (2*se_link)), right_lwr = ilink(fit_link - (2*se_link))) 


ggplot(data=df_meta, aes(x=fit_resp, y=SUM)) +
        geom_point() +
        geom_ribbon(data = df_meta,
                  aes(ymin = right_lwr, ymax = right_upr),
                  alpha = 0.1) 

##################

fit_meta5<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta5, exp=TRUE, digits=4)
summary(fit_meta5)

# Get metrics for model
qchisq(0.95, df.residual(fit_meta5))
deviance(fit_meta5)
pr <- residuals(fit_meta5,"pearson")
sum(pr^2)


fit_meta6<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta5, exp=TRUE, digits=4)

fit_meta7<-glm(SUM~SSU_counts+Age+Gest_age+Antibiotics*Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta7, exp=TRUE, digits=4)

fit_meta8<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta8)
anova(fit_meta8, fit_meta7, test="Chisq")
summ(fit_meta7, exp=TRUE, digits=4)

fit_meta8<-glm(SUM~SSU_counts+Antibiotics+Delivery_mode+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta8, exp=TRUE, digits=4)

fit_meta10<-glm(SUM~Study+SSU_counts+Antibiotics+Delivery_mode+Age+Gest_age+Fortifier+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta10, exp=TRUE, digits=4)

# Tukeys's post hoc test
glht.mod <- glht(fit_meta10, mcp(Study = "Tukey"))
summary(glht(glht.mod))

# Chi squared test
anova(fit_meta8, fit_meta10, test="Chisq")


fit_meta9<-glm(SUM~SSU_counts+Antibiotics+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta9, exp=TRUE, digits=4)

# Chi squared test
anova(fit_meta1, fit_meta2, fit_meta3, fit_meta4,fit_meta9,  fit_meta5, fit_meta7, fit_meta8,test="Chisq")

# Correlation between model variables
cor(model.matrix(fit_meta10)[,-1])

cor(model.matrix(fit_meta6)[,-1])

df_meta$SSU_counts<-df_meta$SSU_counts/1000
df_meta_nonec<-df_meta[!df_meta$Study=="NEC",]
df_meta_nonec$SSU_counts2<-df_meta_nonec$SSU_counts/10000
summ(glmer(SUM~Gest_age+Age+Formula+(1|Study), data=df_meta, family="Gamma"(link="log")), exp=TRUE)


# Test fortifier
fit_meta4_fort<-glm(SUM~Gest_age+Fortifier, data=df_meta_fort, family="Gamma"(link="log") )
summ(fit_meta4_fort, exp=TRUE, digits=4)

# Test in preterm infants (gestatioal age < 37 weeks)
fit_meta4_preterm<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta[df_meta$Gest_age<37,], family="Gamma"(link="log") )
summ(fit_meta4_preterm, exp=TRUE, digits = 4)

fit_meta4_nonec<-glm(SUM~Age+SSU_counts+Gest_age+Formula, data=df_meta[!df_meta$Study%in%c("NEC"),], family="Gamma"(link="log") )
summ(fit_meta4_nonec, exp=TRUE)
summary(fit_meta4_nonec)



fit_meta7<-glm(SUM~SSU_counts+Age+Gest_age*Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta7, exp=TRUE, digits=4)
summary(fit_meta7)

# Chi squared test
anova(fit_meta1, fit_meta2, fit_meta3,fit_meta4, fit_meta5, fit_meta7, test="Chisq")
```

# Plot GLM results
```{r, results="hide", messages =FALSE, warnings =FALSE}
merged_phyloseq_ff<-subset_samples(merged_phyloseq, Study%in%c("Gibson", "NEC"))
merged_phyloseq_ff<-subset_samples(merged_phyloseq_ff, Formula2==FALSE)

df_temp_ff<-data.frame(SUM=sample_sums(merged_phyloseq_ff), Gest_age=sample_data(merged_phyloseq_ff)$Gest_age,
                       Age=sample_data(merged_phyloseq_ff)$Age,
                       Fortifier=sample_data(merged_phyloseq_ff)$Fortifier)

gest.age.metaplot<-ggplot(df_meta, aes(x=Gest_age, y=(SUM), fill=Formula, color=Formula, size=Age)) +
  gamma_smooth()+
  geom_jitter(width = 1, shape=21, color="black")+
    theme_classic()+
    scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes")) +
     scale_color_manual(values=c("darkgoldenrod3", "azure4"), labels=c("No", "Yes") )+
    labs(title="ARGs", x="Gestational age (weeks)", y="ARG sum abundance") +
  theme(axis.title=element_text(size=16), axis.text.x = element_text(size=rel(1.2)), plot.title = element_text(size=16), axis.text=element_text(size=10))+
scale_y_sqrt()+
       guides(fill = guide_legend(override.aes = list(size = 2)),
    size = guide_legend(override.aes = list(linetype = 0)))


gest.age.metaplot

gest.age.plot_ff<-ggplot(df_temp_ff, aes(x=Gest_age, y=SUM, fill=Fortifier, color=Fortifier, size=Age)) +
  gamma_smooth()+
  geom_jitter(width = 1, shape=21, color="black")+
    theme_classic()+
    scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes") ) +
     scale_color_manual(values=c("darkgoldenrod3", "azure4"), labels=c("No", "Yes") )+
    labs(title="ARGs", x="Gestational age (weeks)", y="ARG sum abundance") +
  theme(axis.title=element_text(size=16), axis.text.x = element_text(size=rel(1.2)), plot.title = element_text(size=16), axis.text=element_text(size=10))+
    scale_y_sqrt()+
     guides(fill = guide_legend(override.aes = list(size = 2)),
            size = guide_legend(override.aes = list(linetype = 0)))

```


#Meta-analysis Metaphlan2

Read Metaphlan results from the meta-analysis datasets and find top taxa in each sample.

```{r, results="hide", messages =FALSE, warnings =FALSE}
otu_table_metaphlan<-read.table("merged_abundance_table_species_meta.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
tax_table_metaphlan<-read.table("tax_table_meta.txt", sep="\t", header=FALSE, fill = TRUE, row.names = 1, na.strings = "")
PHY_meta_metaphlan<-phyloseq(otu_table((otu_table_metaphlan), taxa_are_rows=TRUE), sample_data=sample_data(merged_phyloseq), tax_table(as.matrix(tax_table_metaphlan)))

colnames(tax_table(PHY_meta_metaphlan)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

toptaxa_gen<-find.top.taxa(PHY_meta_metaphlan,"Genus")
toptaxa_sp<-find.top.taxa(PHY_meta_metaphlan,"Species")
toptaxa_fam<-find.top.taxa(PHY_meta_metaphlan,"Family")


PHY_meta_metaphlan<-subset_taxa(PHY_meta_metaphlan, Domain=="k__Bacteria")
#Turn into relative data
PHY_meta_metaphlan<-transform_sample_counts(PHY_meta_metaphlan, function(x) x/sum(x))
sample_data(PHY_meta_metaphlan)$Top_gen<-toptaxa_gen$Genus
sample_data(PHY_meta_metaphlan)$Top_sp<-toptaxa_sp$Species
sample_data(PHY_meta_metaphlan)$Top_fam<-toptaxa_fam$Family


PHY_meta_metaphlan<-subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))


sample_data(PHY_meta_metaphlan)$Top_gen[!sample_data(PHY_meta_metaphlan)$Top_gen%in%c("g__Escherichia", "g__Klebsiella", "g__Bifidobacterium", "g__Bacteroides", "g__Staphylococcus", "g__Enterococcus", "g__Enterobacter", "g__Veillonella", "g__Haemophilus", "g__Clostridium")]<-"Other"


# Ordinate

# Square root transform data

temp<-sqrt(otu_table(PHY_meta_metaphlan))

merged_phyloseq_metaphlan_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaphlan)), tax_table(tax_table((PHY_meta_metaphlan))))

# Ordination

PHY_metaphlan_meta_ord<-ordinate(merged_phyloseq_metaphlan_sqrt, method="PCoA", distance="horn")

p<-plot_ordination(merged_phyloseq_metaphlan_sqrt, PHY_metaphlan_meta_ord, color="Top_gen")
p$layers <- p$layers[-1]

metaphlan.meta.ord.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus") + geom_point(size=5, alpha=0.9)  +theme_classic() + labs(title="Microbial community") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, linetype=2, size=1.5,show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=5, color="black") +
   ylim(-0.75, 0.75) +
  xlim(-0.75, 0.75) +
  coord_equal()


#Test differentially abundant species


#Check Enterobacteriaceae differences by explanatory variables with GLMs

df<-data.frame(ENTERO=sample_sums(subset_taxa(PHY_meta_metaphlan, Family=="f__Enterobacteriaceae")), Delivery_mode=sample_data(PHY_meta_metaphlan)$Delivery_mode, AB=sample_data(PHY_meta_metaphlan)$Antibiotics, Fortifier=sample_data(PHY_meta_metaphlan)$Fortifier, Formula=sample_data(PHY_meta_metaphlan)$Formula2, Gest_age=sample_data(PHY_meta_metaphlan)$Gest_age, Age=sample_data(PHY_meta_metaphlan)$Age)

fit<-glm(ENTERO~., data=df, family="quasibinomial")
summ(fit, exp=TRUE, digits=4)

  fit<-glm(ENTERO~Gest_age+Formula, data=df, family="quasibinomial")
  summary(fit)
```

# Metaphlan meta-analysis with DESeq

Get differentially abundant taxa by different diets using DESeq analysis

```{r, results="hide", messages =FALSE, warnings =FALSE}

###DESEQ

#Species and breastfeeding exclusively vs formula
PHY_meta_metaphlan_trim<-prune_taxa(taxa_sums(otu_table(PHY_meta_metaphlan)>0)>nsamples(PHY_meta_metaphlan)/10, PHY_meta_metaphlan)

temp_arg<-otu_table(subset_samples(PHY_meta_metaphlan_trim, !(Formula2==FALSE&Fortifier==TRUE)))

ARG_DSQ<-phyloseq(otu_table(temp_arg, taxa_are_rows=T), sample_data(sample_data(PHY_meta_metaphlan_trim)), tax_table(as.matrix(tax_table(PHY_meta_metaphlan_trim))))

otu_table(ARG_DSQ)<-otu_table(ARG_DSQ)*sample_data(ARG_DSQ)$SSU_counts+1

dds_arg = phyloseq_to_deseq2(ARG_DSQ, ~ Formula2)
dds_arg$Formula2<-as.factor(dds_arg$Formula2)
dds_arg = DESeq(dds_arg, fitType = "mean", test ="Wald")
resultsNames(dds_arg)

res_arg = results(dds_arg, contrast=c("Formula2", "FALSE", "TRUE"))
alpha = 0.05
sigtab_arg = res_arg[which(res_arg$padj < alpha), ]
sigtab_arg = cbind(as(sigtab_arg, "data.frame"), as(tax_table(ARG_DSQ)[rownames(sigtab_arg), ], "matrix"))
sigtab_arg$Formula<-sigtab_arg$log2FoldChange<0
sigtab_arg<-t(as.data.frame(apply(sigtab_arg, 2, function(y) (gsub(".__", "", y)))))
sigtab_arg<-t(as.data.frame(apply(sigtab_arg, 2, function(y) (gsub("_", " ", y)))))
sigtab_arg<-as.data.frame(sigtab_arg)



write.table(sigtab_arg, file="DESEQ_meta_metaphlan.txt", sep="\t", quote=FALSE)

sigtab_arg$log2FoldChange<-as.numeric(as.character(sigtab_arg$log2FoldChange))
reorder(sigtab_arg$Species, sigtab_arg$log2FoldChange)

pal<-wes_palette("IsleofDogs1", 6, type="discrete")
metaphlan.meta.deseq<-ggplot(data = sigtab_arg,
       aes(x = reorder(Species, log2FoldChange),
           y = log2FoldChange)) + 
  geom_bar(stat = "identity", aes(fill = Class), position=position_dodge()) + 
  ylab("Log2 Fold Change") + 
  scale_fill_manual(values=pal)+
  xlab("Species")+
  coord_flip()+
  theme(legend.position = c(0.4,0.2), legend.text = element_text(size=rel(0.7)), axis.text.y = element_text(face="italic") ) +
  guides(fill =guide_legend(ncol=1))
  



```


# Meta-analysis Metaxa2

Metaxa2 result analysis
```{r, results="hide", messages =FALSE, warnings =FALSE}
otu_table_metaxa<-read.table("metaxa_meta-analysis_level6_otutab.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
tax_table_metaxa<-read.table("metaxa_meta-analysis_level6_tax.txt", sep="\t", header=FALSE, fill = TRUE, row.names = 1, na.strings = "")
tax_table_metaxa[is.na(tax_table_metaxa)]<-"Unclassified"

tax_table_metaxa$V8<-paste(tax_table_metaxa$V3, tax_table_metaxa$V4, tax_table_metaxa$V5,  tax_table_metaxa$V6, sep=";")
tax_table_metaxa$V9<-paste(tax_table_metaxa$V3, tax_table_metaxa$V4, tax_table_metaxa$V5,  tax_table_metaxa$V6,  tax_table_metaxa$V7, sep=";")

PHY_meta_metaxa<-phyloseq(otu_table((otu_table_metaxa), taxa_are_rows=TRUE), sample_data=sample_data(merged_phyloseq), tax_table(as.matrix(tax_table_metaxa)))


colnames(tax_table(PHY_meta_metaxa)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Fam_path", "Gen_path")

PHY_meta_metaxa<-subset_taxa(PHY_meta_metaxa, Domain=="Bacteria")
#Turn into relative data
PHY_meta_metaxa_rel<-transform_sample_counts(PHY_meta_metaxa, function(x) x/sum(x))

find.top.taxa <- function(x,taxa){
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- (tax[k,])
  m <- data.frame(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
toptaxa<-find.top.taxa(PHY_meta_metaxa_rel,"Family")

toptaxa_gen<-find.top.taxa(PHY_meta_metaxa_rel,"Genus")


sample_data(PHY_meta_metaxa_rel)$Top_fam_meta<-toptaxa$Family
sample_data(PHY_meta_metaxa_rel)$Top_gen_meta<-paste(toptaxa_gen$Genus, "-" ,toptaxa$Family)
sample_data(PHY_meta_metaxa_rel)$Top_gen_meta[!sample_data(PHY_meta_metaxa_rel)$Top_gen_meta%in%c("Unclassified - Enterobacteriaceae", "Bacteroides - Bacteroidaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Bifidobacterium - Bifidobacteriaceae", "Enterococcus - Enterococcacea", "Clostridium - Clostridiaceae",  "Streptococcus - Streptococcaceae", "Veillonella - Veillonellaceae")]<-"Other"
sample_data(PHY_meta_metaxa_rel)$Top_fam_meta[!sample_data(PHY_meta_metaxa_rel)$Top_fam_meta%in%c("Enterobacteriaceae", "Bacteroidaceae", "Staphylococcaceae", "Bifidobacteriaceae", "Enterococcaceae", "Clostridiaceae")]<-"Other"

```

# Ordination of Metaxa and ARG results
```{r, results="hide", messages =FALSE, warnings =FALSE}
#Ordinate
#Metaxa
temp<-sqrt(otu_table(PHY_meta_metaxa_rel))

merged_phyloseq_metaxa_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaxa_rel)), tax_table(tax_table((PHY_meta_metaxa_rel))))

PHY_SP_meta_ord<-ordinate(merged_phyloseq_metaxa_sqrt, method="PCoA", distance="horn")

p<-plot_ordination(merged_phyloseq_metaxa_sqrt, PHY_SP_meta_ord, color="Top_gen_meta")
p$layers <- p$layers[-1]

metaxa.meta.ord.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", guide=guide_legend(nrow=3),breaks=c("Bacteroides - Bacteroidaceae", "Bifidobacterium - Bifidobacteriaceae", "Clostridium - Clostridiaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Streptococcus - Streptococcaceae", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae", "Other"), labels=c("Bacteroides", "Bifidobacterium", "Clostridium", "Escherichia-Shigella", "Staphylococcus", "Streptococcus", "Uncl. Enterobacteriaceae", "Veillonella", "Other")) + geom_point(size=4)  +theme_classic() + labs(title="Microbial community") + theme(text = element_text(size=16), legend.position = "bottom") +
  stat_ellipse(level = 0.90, linetype=2, size=1.5,show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
   ylim(-0.75, 0.75) +
  xlim(-0.75, 0.75)


#ARGs
temp<-sqrt(otu_table(merged_phyloseq))

merged_phyloseq_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaxa_rel)), tax_table(tax_table((merged_phyloseq))))


PHY_ARG_meta_ord<-ordinate(merged_phyloseq_sqrt, method="PCoA", distance="horn")


p<-plot_ordination(merged_phyloseq_sqrt, PHY_ARG_meta_ord, color="Top_gen_meta")

p$layers <- p$layers[-1]
ARG.meta.ord.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", breaks=c("Bacteroides - Bacteroidaceae", "Bifidobacterium - Bifidobacteriaceae", "Clostridium - Clostridiaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Streptococcus - Streptococcaceae", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae", "Other")) + geom_point(size=4)  +theme_classic() + labs(title="ARGs") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, size=1.5, linetype=2, show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
   ylim(-0.75, 0.75) +
  xlim(-0.75, 0.75) #+
  coord_equal() 


df_temp<-data.frame(SUM=sample_sums(merged_phyloseq), Top_genus=sample_data(merged_phyloseq_sqrt)$Top_gen_meta, Formula=sample_data(merged_phyloseq_sqrt)$Formula2)

fit_meta_genera<-glm(SUM~Top_genus, data=df_temp, family="Gamma"(link="log"))

#Post hoc test for the model
glht.mod <- glht(fit_meta_genera, mcp(Top_genus = "Tukey"))
summary(glht(glht.mod))

# Perform adonis for ARGs

arg_adonis_meta<-pairwise.adonis(x=(t(otu_table(merged_phyloseq_sqrt))), factors = sample_data(merged_phyloseq_sqrt)$Formula2,  sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")
arg_adonis_meta


#Adonis for species
species_adonis_meta<-pairwise.adonis(x=(t(otu_table(merged_phyloseq_metaphlan_sqrt))), factors = sample_data(merged_phyloseq_metaphlan_sqrt)$Formula2,  sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")


```


# Differences in ARG sums related to top genus using GLMs
```{r, results="hide", messages =FALSE, warnings =FALSE}
#Check differences in ARG sums related to top genus using GLMs
df_temp<-data.frame(SUM=sample_sums(merged_phyloseq), Top_genus=sample_data(merged_phyloseq_sqrt)$Top_gen_meta, Formula=sample_data(merged_phyloseq)$Formula2) 


annotation_df <- data.frame(start=c("Bacteroides - Bacteroidaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Other", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae", "Staphylococcus - Staphylococcaceae", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae"), 
                            end=c("Staphylococcus - Staphylococcaceae", "Bifidobacterium - Bifidobacteriaceae", "Bifidobacterium - Bifidobacteriaceae", "Clostridium - Clostridiaceae", "Clostridium - Clostridiaceae", "Escherichia-Shigella - Enterobacteriaceae", "Escherichia-Shigella - Enterobacteriaceae", "Escherichia-Shigella - Enterobacteriaceae", "Other", "Staphylococcus - Staphylococcaceae", "Staphylococcus - Staphylococcaceae"),
                            y=seq(1.1,6.1, by=0.5),
                            label=c("**", "*", "***", "**", "***", "*", "**", "*", "***", "***", "**"))

ARG.gen.sum.plot<-ggplot(df_temp, aes(x=Top_genus, y=log10(SUM)))  +geom_jitter(data=df_temp, aes(x=Top_genus, y=log10((SUM)), color=Top_genus, fill=Top_genus), size=5, width=0.3, shape=21, color="black") + geom_boxplot(alpha=0.5, outlier.shape = NA) +
  scale_fill_brewer("Dominant genus", palette="BrBG") +
   geom_signif(data=annotation_df,
              aes(xmin=start, xmax=end, annotations=label, y_position=y),
              textsize = 4, vjust = 0.3,
              manual=TRUE) +
              ylab("Log10 sum abundance/16S rRNA gene") +
              theme_classic()+
              theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.text = element_text(face="italic")) +
              scale_y_continuous(breaks=seq(-4, 6, 1))
```


# Diversity analysis

```{r, results="hide", messages =FALSE, warnings =FALSE}
# Make dataframe for diversity analysis with GLMs
# Calculate Shannon and Simpson diversity indices
df <- data.frame(DIV=diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "shannon"), DIV2=diversity(t(otu_table(subset_samples(PHY_meta_metaphlan))), index = "simpson"), Formula=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Formula2, Gest_age=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Gest_age,
                 Age=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Age,
                  Antibiotics=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Antibiotics,
                  Delivery_mode=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Delivery_mode)

# Perform GLM analysis
summ(glm(DIV2~Age+Gest_age+Antibiotics+Formula, data=df, family=Gamma(link="log")), exp=TRUE)

df$Formula[df$Formula==TRUE]<-"Yes"
df$Formula[df$Formula==FALSE]<-"No"

# Plot results

metaphlan.shandiv.plot<-ggplot(df, aes(x=Formula, y=DIV, fill=Formula))  +geom_jitter(width=0.4, size=4, shape=21) + geom_boxplot(alpha=0.5) + scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes"))+ theme_classic()  +guides(alpha=FALSE, fill=FALSE, color=FALSE)+
 labs(y="Shannon diversity")  + ggtitle("Shannon diversity of species") 

metaphlan.simdiv.plot<-ggplot(df, aes(x=Formula, y=DIV2, fill=Formula))  +geom_jitter(width=0.4, size=4, shape=21) + geom_boxplot(alpha=0.5) + scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes"))+ theme_classic()  +guides(alpha=FALSE, fill=FALSE, color=FALSE)+
 labs(y="Simpson diversity")  + ggtitle("Simpson diversity of species") 

meta.div.plot<-cowplot::plot_grid(metaphlan.shandiv.plot, metaphlan.simdiv.plot, labels = c("c", "d"))

# Analysis or variance

a0<-aov(diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "shannon")~as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data=as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
TukeyHSD(a0)

a1<-aov(diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "simpson")~as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data=as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
TukeyHSD(a1)


```


# Deseq meta-analysis for ARGS and Metaxa
```{r, results="hide", messages =FALSE, warnings =FALSE}
#ARGs and breastfeeding exclusively vs formula
temp_arg<-otu_table(subset_samples(merged_phyloseq, !(Formula2==FALSE&Fortifier==TRUE)))*5*10^4+1
ARG_tax2<-ARG_tax

# Clean up gene names
ARG_tax2$V3<-gsub('TEM-[[:digit:]]+', 'TEM', ARG_tax2$V3)
ARG_tax2$V3<-gsub('TEMA', 'TEM', ARG_tax2$V3)
ARG_tax2$V3<-gsub('TEMB', 'TEM', ARG_tax2$V3)
ARG_tax2$V3<-gsub('TEMC', 'TEM', ARG_tax2$V3)
ARG_tax2$V3<-gsub('erm\\(C\\)', 'ermC', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ErmC', 'ermC', ARG_tax2$V3)
ARG_tax2$V3<-gsub('erm\\(B\\)', 'ermB', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ErmB', 'ermB', ARG_tax2$V3)
ARG_tax2$V3<-gsub('Cfx', 'cfx', ARG_tax2$V3)
ARG_tax2$V3<-gsub('APH', 'aph', ARG_tax2$V3)
ARG_tax2$V3<-gsub('FosA', 'fosA', ARG_tax2$V3)
ARG_tax2$V3<-gsub('AAC', 'aac', ARG_tax2$V3)
ARG_tax2$V3<-gsub('ANT', 'ant', ARG_tax2$V3)

# DESeq analysis
ARG_DSQ<-phyloseq(otu_table(temp_arg, taxa_are_rows=T), sample_data(sample_data(merged_phyloseq)), tax_table(as.matrix(ARG_tax2)))
ARG_DSQ_glom<-tax_glom(ARG_DSQ, taxrank="V3")
dds_arg = phyloseq_to_deseq2(ARG_DSQ_glom, ~ Formula2)
dds_arg = DESeq(dds_arg, fitType = "mean", test ="Wald")
resultsNames(dds_arg)

res_arg = results(dds_arg)
alpha = 0.05
sigtab_arg = res_arg[which(res_arg$padj < alpha), ]
sigtab_arg = cbind(as(sigtab_arg, "data.frame"), as(tax_table(ARG_DSQ)[rownames(sigtab_arg), ], "matrix"))
sigtab_arg<-sigtab_arg[sigtab_arg$baseMean>5,]

otu_table(ARG_DSQ_glom)[otu_table(ARG_DSQ_glom)==1]<-0
otu_table(ARG_DSQ_glom)[otu_table(ARG_DSQ_glom)>0]<-1
n<-rowSums(otu_table(ARG_DSQ))

sigtab_arg=merge(sigtab_arg, as.data.frame(n), by=0)

sigtab_arg<-sigtab_arg[sigtab_arg$n>20,]

# Plot results
Deseq_meta<-ggplot(sigtab_arg, aes(x=V3, y=log2FoldChange, color=V2, fill=V2, size=100*baseMean)) + geom_jitter(alpha=0.8, width = 0.1, shape=21, color="black")+theme_minimal()+ theme(axis.title.y = element_text(size=rel(1.3)), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)), panel.grid.major = element_blank(), legend.position="bottom", axis.text.x=element_blank())+ geom_text_repel(aes(label=V3), size=4, color= "black", fontface = "italic") +  scale_fill_brewer(palette="BrBG", "ARG class") + scale_size(range=c(3,20),breaks=c(3, 5, 10, 15, 20)) + labs(color="") +guides(size=FALSE) +  ylim(-6, 6) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") + xlab("") +
guides(fill = guide_legend(override.aes = list(size = 4)),
            size = guide_legend(override.aes = list(linetype = 0))) +
  guides(size = guide_legend(override.aes = list(linetype = 0)))

# Metaxa DESEq analysis

# Taxa and breastfeeding exclusively vs formula
temp_mtx<-otu_table(PHY_meta_metaxa)+1
mtx_tax2<-tax_table(PHY_meta_metaxa)

mtx_DSQ<-phyloseq(otu_table(temp_mtx, taxa_are_rows=T), sample_data(PHY_meta_metaxa), tax_table(as.matrix(mtx_tax2)))
mtx_DSQ_glom<-tax_glom(mtx_DSQ, taxrank="Genus")
dds_mtx = phyloseq_to_deseq2(mtx_DSQ_glom, ~ Formula2)
dds_mtx = DESeq(dds_mtx, fitType = "mean", test ="Wald")
res_mtx = results(dds_mtx, cooksCutoff = FALSE)
alpha = 0.05
sigtab_mtx = res_mtx[which(res_mtx$padj < alpha), ]
sigtab_mtx = cbind(as(sigtab_mtx, "data.frame"), as(tax_table(mtx_DSQ)[rownames(sigtab_mtx), ], "matrix"))

sigtab_mtx<-sigtab_mtx[1:12]
sigtab_mtx<-sigtab_mtx[sigtab_mtx$baseMean>10,]
sigtab_mtx<-sigtab_mtx[!sigtab_mtx$Family%in%c("Unclassified", "Family"),]
kable(sigtab_mtx, caption="Genera different by formula")


# Plot results

pal<-wes_palette("IsleofDogs1", 6, type="discrete")

Deseq_meta_mtx<-ggplot(sigtab_mtx, aes(x=Family, y=log2FoldChange, color=Class, fill=Class, size=100*baseMean)) + geom_jitter(alpha=0.9, width = 0.2, shape=21, color="black")+theme_minimal()+ theme(axis.text.x = element_text(angle =45, hjust = 0.9, vjust = 0.97, size=rel(1.3), face="italic"), axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_text(size=rel(1.3)), legend.text = element_text(size=rel(1.2)), legend.title = element_text(size=rel(1.2)), legend.position="bottom")+ geom_text_repel(aes(label=Genus, fontface="italic"), size=4, color= "black") +  scale_fill_manual(values = pal) + scale_size(range=c(3,20),breaks=c(3, 5, 10, 15, 20)) + labs(color="") +guides(size=FALSE) +  ylim(-5, 5) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") +
guides(fill = guide_legend(override.aes = list(size = 4)),
            size = guide_legend(override.aes = list(linetype = 0)))


```

# Nature data set analyses
https://doi.org/10.1038/s41586-019-1560-1

```{r, results="hide", messages =FALSE, warnings =FALSE}
Nature<-read.csv("Nature_metadata.csv", header=TRUE, sep=";")
Nature$Individual<-gsub( "_T[[:digit:]]+", "", Nature$Individual)
table(Nature$Time_point)
unique(Nature$Individual)
table(Nature$Feeding_method)
Nature$Formula<-!Nature$Feeding_method=="BF"

Nature_uniq<-Nature[(match(unique(Nature$Individual), Nature$Individual)),]
write.table(Nature_uniq, file="Nature_fist_samples.csv", sep = ";", quote = FALSE)

Nature_metadata<-read.table(file="Nature_fist_samples.csv", sep = ";", header = TRUE, row.names = 2)

# ARGs

Nature_ARG_bt_res<-as.matrix(read.table("ARG_genemat_Nature.txt", fill= 1, header= T, row.names = 1, check.names = T))
ARG_tax<-read.table("ARG_tax_table_final.txt", fill=1, row.names = 1, header=F)
ARG_tax$V2<-gsub('_resistance', '', ARG_tax$V2)
ARG_tax$V2<-gsub('_', ' ', ARG_tax$V2)
ARG_lengths_res <- as.matrix(read.table("ARG_genelenghts_resfinder.txt", fill= 1, header= T, row.names = 1, check.names = F))
ARG_lengths<- as.matrix(read.table("ARG_genelenghts.txt", fill= 1, header= T, row.names = 1, check.names = F))

# Read in 16S counts

#Read in SSU_counts
Nature_SSU_counts <- as.matrix(read.table("Nature_SSU_counts", fill= 1, header= F, row.names = 1, check.names = F))

#Divide by genelengths
#Divide by ARG gene lengths
Nature_ARG_bt_res<-Nature_ARG_bt_res[row.names(Nature_ARG_bt_res)%in%row.names(ARG_lengths_res),]

Nature_ARG_bt_res<-Nature_ARG_bt_res[,colnames(Nature_ARG_bt_res)%in%c(rownames(Nature_SSU_counts))]

Nature_arg_length_norm_res <- Nature_ARG_bt_res/ARG_lengths_res[,1]

#Read in SSU counts and tax table
Nature_ARG_SSU_length_norm_res<-t(t(Nature_arg_length_norm_res)/Nature_SSU_counts[,1])*1541

Nature_ARG_SSU_length_norm_res[is.na(Nature_ARG_SSU_length_norm_res)]<-0

#Make phyloseq
Nature_ARG_PHY_res<-phyloseq(otu_table(Nature_ARG_SSU_length_norm_res, taxa_are_rows = T), sample_data(Nature_metadata), tax_table(as.matrix(ARG_tax)))

sample_data(Nature_ARG_PHY_res)$SSU_counts<-Nature_SSU_counts[,1]

Nature_ARG_PHY_res<-subset_samples(Nature_ARG_PHY_res, !colSums(otu_table(Nature_ARG_PHY_res))=="NA"&colSums(otu_table(Nature_ARG_PHY_res))>0)

Nature_ARG_PHY_res<-subset_samples(Nature_ARG_PHY_res, SSU_counts>1000)

df_Nature<-data.frame(SUM=(sample_sums(Nature_ARG_PHY_res)), 
                Formula=as.factor(sample_data(Nature_ARG_PHY_res)$Formula),
                 Age=sample_data(Nature_ARG_PHY_res)$Time_point,
                Delivery_mode=sample_data(Nature_ARG_PHY_res)$Delivery_mode,
                 Country="UK",
                Study="Nature",
                Gender=sample_data(Nature_ARG_PHY_res)$Gender,
                SSU_counts=sample_data(Nature_ARG_PHY_res)$SSU_counts,
                Antibiotics=sample_data(Nature_ARG_PHY_res)$Abx_Baby_in_hospital,
                Antibiotics2=sample_data(Nature_ARG_PHY_res)$Abx_Baby_after_hospital,
                Antibiotics_mom=sample_data(Nature_ARG_PHY_res)$Abx_mother_prior_birth,
                IAP=sample_data(Nature_ARG_PHY_res)$Abx_mother_labour_IAP,
                Weight=sample_data(Nature_ARG_PHY_res)$Birth_weight,
                Diet=sample_data(Nature_ARG_PHY_res)$Feeding_method,
                Breast_feeding=sample_data(Nature_ARG_PHY_res)$Breastfeeding_status,
                Breast_feeding_1hr=sample_data(Nature_ARG_PHY_res)$Breastfeeding_1hr_birth,
                Hospital=sample_data(Nature_ARG_PHY_res)$Hospital)

sample_data(Nature_ARG_PHY_res)$Age<-sample_data(Nature_ARG_PHY_res)$Time_point
sample_data(Nature_ARG_PHY_res)$Antibiotics<-sample_data(Nature_ARG_PHY_res)$Abx_Baby_in_hospital=="Yes"
sample_data(Nature_ARG_PHY_res)$SUM<-sample_sums(Nature_ARG_PHY_res)
sample_data(Nature_ARG_PHY_res)$Formula2<-sample_data(Nature_ARG_PHY_res)$Formula
sample_data(Nature_ARG_PHY_res)$Delivery_mode<-gsub("Caesarean", "CS", sample_data(Nature_ARG_PHY_res)$Delivery_mode) 
sample_data(Nature_ARG_PHY_res)$Delivery_mode<-gsub("Vaginal", "V", sample_data(Nature_ARG_PHY_res)$Delivery_mode)
sample_data(Nature_ARG_PHY_res)$Study<-"Nature"
sample_data(Nature_ARG_PHY_res)$Country<-"UK"
sample_data(Nature_ARG_PHY_res)$Donor_milk<-0
sample_data(Nature_ARG_PHY_res)$Fortifier<-"FALSE"
sample_data(Nature_ARG_PHY_res)$Gest_age<-39.5

# Metaphlan

otu_table_metaphlan_nature<-read.table("merged_abundance_table_species_Nature.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
tax_table_metaphlan_nature<-read.table("tax_table_species_Nature.txt", sep="\t", header=FALSE, fill = TRUE, row.names = 1, na.strings = "")
PHY_nature_metaphlan<-phyloseq(otu_table((otu_table_metaphlan_nature), taxa_are_rows=TRUE), sample_data=sample_data(Nature_ARG_PHY_res), tax_table(as.matrix(tax_table_metaphlan_nature)))

colnames(tax_table(PHY_nature_metaphlan)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

toptaxa_gen<-find.top.taxa(PHY_nature_metaphlan,"Genus")
toptaxa_sp<-find.top.taxa(PHY_nature_metaphlan,"Species")
toptaxa_fam<-find.top.taxa(PHY_nature_metaphlan,"Family")

PHY_nature_metaphlan<-subset_taxa(PHY_nature_metaphlan, Domain=="k__Bacteria")
# Turn into relative data
PHY_nature_metaphlan<-transform_sample_counts(PHY_nature_metaphlan, function(x) x/sum(x))
sample_data(PHY_nature_metaphlan)$Top_gen<-toptaxa_gen$Genus
sample_data(PHY_nature_metaphlan)$Top_sp<-toptaxa_sp$Species
sample_data(PHY_nature_metaphlan)$Top_fam<-toptaxa_fam$Family

Nature_ARG_PHY_res_metaphlan<-subset_samples(Nature_ARG_PHY_res, sample_names(Nature_ARG_PHY_res)%in%c(sample_names(PHY_nature_metaphlan)))
PHY_nature_metaphlan_res<-subset_samples(PHY_nature_metaphlan, sample_names(PHY_nature_metaphlan)%in%c(sample_names(Nature_ARG_PHY_res)))

sample_data(Nature_ARG_PHY_res_metaphlan)$Top_gen<-sample_data(PHY_nature_metaphlan_res)$Top_gen
df_temp<-data.frame(SUM=sample_sums(Nature_ARG_PHY_res_metaphlan),Top_genus=sample_data(Nature_ARG_PHY_res_metaphlan)$Top_gen)

sample_data(Nature_ARG_PHY_res_metaphlan)$Top_gen[!sample_data(Nature_ARG_PHY_res_metaphlan)$Top_gen%in%c("g__Escherichia", "g__Klebsiella", "g__Bifidobacterium", "g__Bacteroides", "g__Staphylococcus", "g__Enterococcus", "g__Enterobacter", "g__Veillonella", "g__Haemophilus", "g__Clostridium", "g__Streptococcus")]<-"Other"

df_temp<-data.frame(SUM=sample_sums(Nature_ARG_PHY_res_metaphlan), Top_genus=sample_data(Nature_ARG_PHY_res_metaphlan)$Top_gen)

ggplot(df_temp, aes(x=Top_genus, y=log10(SUM)))  +geom_jitter(data=df_temp, aes(x=Top_genus, y=log10((SUM)), color=Top_genus, fill=Top_genus), size=5, width=0.3, shape=21, color="black") + geom_boxplot(alpha=0.5, outlier.shape = NA) 


#Merge phyloseq
merged_phyloseq_Nature<-merge_phyloseq(merged_phyloseq, Nature_ARG_PHY_res_metaphlan)


sample_data(merged_phyloseq_Nature)<-sample_data(merged_phyloseq_Nature)[,colnames(sample_data(merged_phyloseq_Nature))%in%c("Formula2", "Delivery_mode", "Age", "Antibiotics", "Gest_age", "SSU_counts", "Study", "Sample", "Country", "Donor_milk", "Twin", "Fortifier")]


# Make dataframe

df_meta_Nature<-data.frame(SUM=sample_sums(merged_phyloseq_Nature),
                    Formula=as.factor(sample_data(merged_phyloseq_Nature)$Formula2),
                    Fortifier=(sample_data(merged_phyloseq_Nature)$Fortifier),
                    Antibiotics=sample_data(merged_phyloseq_Nature)$Antibiotics, Delivery_mode=sample_data(merged_phyloseq_Nature)$Delivery_mode,
                    Age=sample_data(merged_phyloseq_Nature)$Age,
                    Gest_age=sample_data(merged_phyloseq_Nature)$Gest_age,
                    Study=sample_data(merged_phyloseq_Nature)$Study,
                    Sample=sample_names(merged_phyloseq_Nature),
                    Country=sample_data(merged_phyloseq_Nature)$Country,
                    SSU_counts=sample_data(merged_phyloseq_Nature)$SSU_counts,
                    Twin=sample_data(merged_phyloseq_Nature)$Twin)


df_meta_Nature_NFF<-df_meta_Nature[df_meta_Nature$Fortifier=="FALSE",]

df_meta_Nature_F<-df_meta_Nature[df_meta_Nature$Formula==TRUE,]
df_meta_Nature_BM<-df_meta_Nature[df_meta_Nature$Formula==FALSE,]

table(df_meta_Nature_F$Antibiotics)
table(df_meta_Nature_BM$Antibiotics)
table(df_meta_Nature_F$Delivery_mode)
table(df_meta_Nature_BM$Delivery_mode)
table(df_meta_Nature_BM$Country)
table(df_meta_Nature_F$Country)
mean(df_meta_Nature_BM$Age)
mean(df_meta_Nature_F$Age)
mean(df_meta_Nature_F$Gest_age)
mean(df_meta_Nature_BM$Gest_age)

fit0<-(glm(SUM~Country+Age+Gest_age+Delivery_mode+Antibiotics+Gest_age+Formula, data=df_meta_Nature, family="Gamma"(link="log")))
summ(fit0, exp=TRUE)

fit1<-(glm(SUM~Country+Age+Gest_age+Delivery_mode+Antibiotics+Gest_age*Formula, data=df_meta_Nature, family="Gamma"(link="log")))
summ(fit1, exp=TRUE)

fit2<-(glm(SUM~Study+Age+Gest_age+Delivery_mode+Antibiotics+Gest_age*Formula, data=df_meta_Nature, family="Gamma"(link="log")))
summ(fit2, exp = TRUE)

anova(fit0, fit1, fit2, test="Chisq")

cor(model.matrix(fit1)[,-1])

modelEffectSizes(fit1)
modelPower(pc=5, pa=6, peta2 = 0.006,power=0.5)
```

# Function for testing similarity between twins and mother-infant pairs
```{r, results="hide", messages =FALSE, warnings =FALSE}
#Read in function for testing significances
lmp = function(y,x,n.perms = 9999){
   test = numeric(n.perms)
   test[1] = abs(summary(lm(y~x))$coefficients[2,3])
   for(ii in 2:n.perms){
       test[ii] = abs(summary(lm(sample(y)~x))$coefficients[2,3])
   }
   p.value = sum(test>=test[1])/n.perms
   return(p.value)
}
```

#Session Info

```{r}
sessionInfo()
```

