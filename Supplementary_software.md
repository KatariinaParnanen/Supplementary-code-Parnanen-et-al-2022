Supplementary Code
================
Katariina Pärnänen
Sep/13/2019

## Set environment

Acquire required packages

``` r
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
##
```

## Read in metaphlan results on species level

The abundance\_species.txt file is created from the merged abundance
table from Metaphlan. Lines with ’s\_\_’ are picked and the lines with
’t\_\_’ are removed. The abundance table then has only species level
entries. Column with OTU and running number is added after
this.

``` r
metaphlan_sp <- as.matrix(read.table("merged_abundance_table_species.txt", fill= 1, header= T, row.names = 1, check.names = F))
```

## Read in taxonomy table for species level

Taxonomy table is created from the species level merged abundance table
with awk and sed scripts.

``` r
tax_sp<- read.table(("tax_table_species.txt"), fill=1, row.names=1)
tax_sp<-apply(tax_sp, 2, function(y) (gsub(".__", "", y)))
```

## Read in metadata

``` r
sample_data<-read.csv(as.matrix("ForAnalysis_NEC_metadata_46.csv"), header=T, row.names = 1,sep = ";", stringsAsFactors=FALSE)
sample_data$Pair<-substring(row.names(sample_data), 1, 6)
sample_data$Sample<-row.names(sample_data)
sample_data$Gest_age2<-cut(sample_data$Gest_age, breaks=3)
sample_data$Antibiotics<-sample_data$Inf_AB=="Y"
sample_data$Donor_Milk<-sample_data$Donor_Breast_Milk=="y"
```

## Merge into a phyloseq object

``` r
PHY_SP <- phyloseq(otu_table(metaphlan_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_sp)), sample_data(sample_data))
```

## Change the taxonomic levels to something meaningful

``` r
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

``` r
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

``` r
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
```

    ## [1] TRUE

``` r
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

``` r
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

``` r
human_sp <- as.matrix(read.table("all_genefamilies_cpmEC_genemat.tsv", fill= 1, header= T, check.names = F, sep = "\t"))
OTU<- paste("OTU", 1:nrow(human_sp), sep="")
rownames(human_sp)<- OTU
```

## Read in taxonomy table of genes

``` r
tax_EC<- read.table(("EC_taxtable.txt"), fill=TRUE, header=TRUE, sep ="\t", row.names = 1, quote =  "")
OTU<- paste("OTU", 1:nrow(tax_EC ), sep="")
rownames(tax_EC )<- OTU
```

## Merge into a phyloseq object

``` r
PHY_humann <- phyloseq(otu_table(human_sp, taxa_are_rows = TRUE), tax_table(as.matrix(tax_EC)), sample_data(sample_data))
```

## Ordinate and calculate diversities

``` r
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
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

    ##                        pairs   F.Model        R2 p.value p.adjusted sig
    ## 1        Klebsiella vs Other  5.001299 0.1562843   0.001     0.0012   *
    ## 2  Klebsiella vs Escherichia 19.529024 0.4591929   0.001     0.0012   *
    ## 3  Klebsiella vs Veillonella  8.256940 0.2922093   0.001     0.0012   *
    ## 4       Other vs Escherichia  6.588992 0.2304730   0.001     0.0012   *
    ## 5       Other vs Veillonella  3.006875 0.1366334   0.005     0.0050   *
    ## 6 Escherichia vs Veillonella 19.957733 0.5709104   0.001     0.0012   *

``` r
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
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

    ##                        pairs  F.Model         R2 p.value p.adjusted sig
    ## 1        Klebsiella vs Other 6.589081 0.19616735   0.001     0.0030   *
    ## 2  Klebsiella vs Escherichia 8.589191 0.27190285   0.001     0.0030   *
    ## 3  Klebsiella vs Veillonella 4.428655 0.18128935   0.002     0.0040   *
    ## 4       Other vs Escherichia 1.646306 0.06962212   0.043     0.0645    
    ## 5       Other vs Veillonella 1.231788 0.06088381   0.211     0.2110    
    ## 6 Escherichia vs Veillonella 1.337707 0.08187849   0.165     0.1980

``` r
# Check separation in MGEs

pairwise.adonis(x=sqrt(t(otu_table(MGE_PHY))), factors=sample_data(MGE_PHY)$Top_genus2, sim.function='vegdist', sim.method='horn',p.adjust.m='fdr')
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

    ##                        pairs  F.Model         R2 p.value p.adjusted sig
    ## 1        Klebsiella vs Other 1.945615 0.06721624   0.023     0.0460   .
    ## 2  Klebsiella vs Escherichia 4.530901 0.16457510   0.001     0.0060   *
    ## 3  Klebsiella vs Veillonella 1.615640 0.07474403   0.071     0.0852    
    ## 4       Other vs Escherichia 2.353726 0.09664748   0.007     0.0210   .
    ## 5       Other vs Veillonella 1.329625 0.06540330   0.154     0.1540    
    ## 6 Escherichia vs Veillonella 2.400935 0.13797733   0.031     0.0465   .

``` r
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
```

    ## <ggproto object: Class CoordFixed, CoordCartesian, Coord, gg>
    ##     aspect: function
    ##     backtransform_range: function
    ##     clip: on
    ##     default: FALSE
    ##     distance: function
    ##     expand: TRUE
    ##     is_free: function
    ##     is_linear: function
    ##     labels: function
    ##     limits: list
    ##     modify_scales: function
    ##     range: function
    ##     ratio: 1
    ##     render_axis_h: function
    ##     render_axis_v: function
    ##     render_bg: function
    ##     render_fg: function
    ##     setup_data: function
    ##     setup_layout: function
    ##     setup_panel_params: function
    ##     setup_params: function
    ##     transform: function
    ##     super:  <ggproto object: Class CoordFixed, CoordCartesian, Coord, gg>

``` r
# Check significance of Humann clusters
pairwise.adonis(sqrt(t(otu_table(PHY_humann_gene_prop))), factors = sample_data(PHY_humann_gene_prop)$Top_genus2, sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

    ##                        pairs   F.Model         R2 p.value p.adjusted sig
    ## 1        Klebsiella vs Other  5.318223 0.16981240   0.001     0.0015   *
    ## 2  Klebsiella vs Escherichia 25.987757 0.54154973   0.001     0.0015   *
    ## 3  Klebsiella vs Veillonella  7.341633 0.27870835   0.001     0.0015   *
    ## 4       Other vs Escherichia  9.368032 0.29864903   0.001     0.0015   *
    ## 5       Other vs Veillonella  2.077923 0.09858289   0.032     0.0320   .
    ## 6 Escherichia vs Veillonella 12.884016 0.46205740   0.002     0.0024   *

``` r
# Diversity
df <- data.frame(DIV=diversity(t(otu_table(PHY_SP)), index = "shannon"),
                 DIV2=colSums(otu_table(PHY_SP)>0),
                  Top_genus2=sample_data(PHY_SP)$Top_genus2, Top_genus2=sample_data(PHY_SP)$Top_genus2)

metaphlan.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(alpha=FALSE)+
theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("Microbiota")

# Shannon diversity 
a0<-aov(diversity(t(otu_table(PHY_SP)), index = "shannon")~sample_data(PHY_SP)$Top_genus2, data=as.data.frame(otu_table(PHY_SP)))
TukeyHSD(a0)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity(t(otu_table(PHY_SP)), index = "shannon") ~ sample_data(PHY_SP)$Top_genus2, data = as.data.frame(otu_table(PHY_SP)))
    ## 
    ## $`sample_data(PHY_SP)$Top_genus2`
    ##                                diff        lwr       upr     p adj
    ## Klebsiella-Escherichia  -0.02413604 -0.6417302 0.5934581 0.9995842
    ## Other-Escherichia        0.14549150 -0.4808633 0.7718463 0.9246709
    ## Veillonella-Escherichia  0.10303076 -0.6424806 0.8485421 0.9825114
    ## Other-Klebsiella         0.16962754 -0.3925431 0.7317982 0.8507143
    ## Veillonella-Klebsiella   0.12716680 -0.5652945 0.8196281 0.9606213
    ## Veillonella-Other       -0.04246074 -0.7427467 0.6578252 0.9984567

``` r
df <- data.frame(DIV=diversity(t(otu_table(PHY_humann_gene_prop)), index = "shannon"), DIV2=colSums(otu_table(PHY_humann_gene_prop)>0), Top_genus2=sample_data(PHY_humann_gene_prop)$Top_genus2, Top_genus2=sample_data(PHY_humann_gene_prop)$Top_genus2)
humann.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("Functional genes")

a0<-aov(diversity(t(otu_table(PHY_humann_gene_prop)), index = "shannon")~sample_data(PHY_humann_gene_prop)$Top_genus2, data=as.data.frame(otu_table(PHY_humann_gene_prop)))
TukeyHSD(a0)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity(t(otu_table(PHY_humann_gene_prop)), index = "shannon") ~ sample_data(PHY_humann_gene_prop)$Top_genus2, data = as.data.frame(otu_table(PHY_humann_gene_prop)))
    ## 
    ## $`sample_data(PHY_humann_gene_prop)$Top_genus2`
    ##                                diff         lwr       upr     p adj
    ## Klebsiella-Escherichia  -0.03797296 -0.37452035 0.2985744 0.9902775
    ## Other-Escherichia        0.09485743 -0.24168996 0.4314048 0.8741298
    ## Veillonella-Escherichia  0.48627628  0.08570477 0.8868478 0.0118789
    ## Other-Klebsiella         0.13283039 -0.17439394 0.4400547 0.6564404
    ## Veillonella-Klebsiella   0.52424924  0.14797782 0.9005207 0.0031303
    ## Veillonella-Other        0.39141885  0.01514743 0.7676903 0.0387326

``` r
# Calculate diverisity and make plot

df <- data.frame(DIV=diversity(t(otu_table(ARG_PHY_res2)), index = "shannon"),
                 DIV2=colSums(otu_table(ARG_PHY_res2)>0),
                 Top_genus2=sample_data(ARG_PHY_res2)$Top_genus2, Top_genus2=sample_data(ARG_PHY_res2)$Top_genus2)

ARG.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG", "Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank())  +labs(y="Shannon diversity", x='') + ggtitle("ARGs")

# Analysis of variance

a0<-aov(diversity((t(otu_table(ARG_PHY_res2))), index = "shannon")~sample_data(ARG_PHY_res2)$Top_genus2, data=as.data.frame(otu_table(ARG_PHY_res2)))
TukeyHSD(a0)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity((t(otu_table(ARG_PHY_res2))), index = "shannon") ~ sample_data(ARG_PHY_res2)$Top_genus2, data = as.data.frame(otu_table(ARG_PHY_res2)))
    ## 
    ## $`sample_data(ARG_PHY_res2)$Top_genus2`
    ##                               diff        lwr       upr     p adj
    ## Klebsiella-Escherichia  -0.3034547 -1.0632780 0.4563685 0.7104430
    ## Other-Escherichia       -0.5546926 -1.3252940 0.2159089 0.2330863
    ## Veillonella-Escherichia -0.7520183 -1.6692175 0.1651809 0.1417753
    ## Other-Klebsiella        -0.2512379 -0.9428738 0.4403981 0.7661592
    ## Veillonella-Klebsiella  -0.4485636 -1.3004955 0.4033684 0.5012740
    ## Veillonella-Other       -0.1973257 -1.0588843 0.6642329 0.9274891

``` r
# Calculate diverisity and make plot

df <- data.frame(DIV=diversity(t(otu_table(MGE_PHY)), index = "shannon"),
                 DIV2=colSums(otu_table(MGE_PHY)>0),
                 Top_genus2=sample_data(MGE_PHY)$Top_genus2, Top_genus2=sample_data(MGE_PHY)$Top_genus2)

MGE.div.plot<-ggplot(df, aes(x=Top_genus2, y=DIV, fill=Top_genus2, alpha=0.5)) + geom_boxplot()  + scale_fill_brewer(palette = "BrBG","Dominant genus")+ theme_classic() +guides(fill=FALSE, alpha=FALSE)+
  theme(axis.text.x=element_blank()) +labs(y="Shannon diversity", x='') + ggtitle("MGEs")

# Analysis of variance
a0<-aov(diversity((t(otu_table(MGE_PHY))), index = "shannon")~sample_data(MGE_PHY)$Top_genus2, data=as.data.frame(otu_table(MGE_PHY)))
TukeyHSD(a0)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity((t(otu_table(MGE_PHY))), index = "shannon") ~ sample_data(MGE_PHY)$Top_genus2, data = as.data.frame(otu_table(MGE_PHY)))
    ## 
    ## $`sample_data(MGE_PHY)$Top_genus2`
    ##                                diff        lwr       upr     p adj
    ## Klebsiella-Escherichia   0.44476343 -0.5249177 1.4144445 0.6135470
    ## Other-Escherichia       -0.20189060 -1.1853267 0.7815455 0.9462685
    ## Veillonella-Escherichia  0.42264624 -0.7478769 1.5931694 0.7693962
    ## Other-Klebsiella        -0.64665403 -1.5293150 0.2360069 0.2195138
    ## Veillonella-Klebsiella  -0.02211719 -1.1093467 1.0651123 0.9999412
    ## Veillonella-Other        0.62453684 -0.4749782 1.7240519 0.4352641

``` r
# Cowplot

legend<-get_legend(plot=metaphlan.plot)

cowplot::plot_grid(metaphlan.plot +theme(legend.position = "none"),
                           humann.plot+theme(legend.position = "none"),
                           ARG.plot+theme(legend.position = "none"), 
                           MGE.plot+theme(legend.position = "none"),
                           labels = c("auto"), align = "hv", ncol = 2)
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Random forests for data exploration to pick best predictors for sum abundances

``` r
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

    ## rf variable importance
    ## 
    ##              Overall
    ## Gest_age      100.00
    ## Diet           96.59
    ## M.Seqs         48.12
    ## Age_sampling   40.15
    ## Inf_AB         38.91
    ## Mat_AB         20.73
    ## Inf_inf        19.94
    ## Delivery_met    0.00

## Build model for ARG and MGE sum abundance using gamma distributed GLMs

``` r
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
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ Gest_age + DIET_CODES_4, family = poisson, 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -104.228   -63.141   -25.106     7.196   278.450  
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    12.694053   0.042380  299.53   <2e-16 ***
    ## Gest_age       -0.159377   0.001302 -122.45   <2e-16 ***
    ## DIET_CODES_4F   1.556552   0.011015  141.31   <2e-16 ***
    ## DIET_CODES_4FF  0.788821   0.010912   72.29   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 322757  on 45  degrees of freedom
    ## Residual deviance: 282824  on 42  degrees of freedom
    ## AIC: 283286
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
qchisq(0.95, df.residual(mp))
```

    ## [1] 58.12404

``` r
deviance(mp)
```

    ## [1] 282824.2

``` r
pr <- residuals(mp,"pearson")

sum(pr^2)
```

    ## [1] 392391.6

``` r
#95% cutoff value for chisq
phi <- sum(pr^2)/df.residual(mp)


round(c(phi,sqrt(phi)),4)
```

    ## [1] 9342.6572   96.6574

``` r
# Quasipoisson

mq <- glm(SUM~Gest_age+DIET_CODES_4, family=quasipoisson, data=df)

summary(mq)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ Gest_age + DIET_CODES_4, family = quasipoisson, 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -104.228   -63.141   -25.106     7.196   278.450  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)     12.6941     4.0964   3.099  0.00346 **
    ## Gest_age        -0.1594     0.1258  -1.267  0.21219   
    ## DIET_CODES_4F    1.5566     1.0647   1.462  0.15120   
    ## DIET_CODES_4FF   0.7888     1.0548   0.748  0.45871   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 9342.657)
    ## 
    ##     Null deviance: 322757  on 45  degrees of freedom
    ## Residual deviance: 282824  on 42  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
se <- function(model) sqrt(diag(vcov(model)))

round(data.frame(p=coef(mp), q=coef(mq), se.p=se(mp), se.q=se(mq), ratio=se(mq)/se(mp)), 4)
```

    ##                      p       q   se.p   se.q   ratio
    ## (Intercept)    12.6941 12.6941 0.0424 4.0964 96.6574
    ## Gest_age       -0.1594 -0.1594 0.0013 0.1258 96.6574
    ## DIET_CODES_4F   1.5566  1.5566 0.0110 1.0647 96.6574
    ## DIET_CODES_4FF  0.7888  0.7888 0.0109 1.0548 96.6574

``` r
# Negative binomial

mnb <- glm.nb(SUM~Gest_age+DIET_CODES_4, data=df)

summary(mnb)
```

    ## 
    ## Call:
    ## glm.nb(formula = SUM ~ Gest_age + DIET_CODES_4, data = df, init.theta = 1.085231823, 
    ##     link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.1029  -1.2086  -0.3955   0.1992   2.4959  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    15.94164    2.85579   5.582 2.37e-08 ***
    ## Gest_age       -0.25515    0.08504  -3.000   0.0027 ** 
    ## DIET_CODES_4F   1.49787    0.47995   3.121   0.0018 ** 
    ## DIET_CODES_4FF  0.50606    0.51636   0.980   0.3271    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(1.0852) family taken to be 1)
    ## 
    ##     Null deviance: 63.663  on 45  degrees of freedom
    ## Residual deviance: 52.598  on 42  degrees of freedom
    ## AIC: 889.42
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  1.085 
    ##           Std. Err.:  0.201 
    ## 
    ##  2 x log-likelihood:  -879.421

``` r
# Variance
1/mnb$theta
```

    ## [1] 0.9214621

``` r
mnbg <- glm(SUM~Gest_age+DIET_CODES_4, family=negative.binomial(mnb$theta), data=df)

all(abs(coef(mnb)==coef(mnbg)) < 1e-12)
```

    ## [1] TRUE

``` r
v = 1/mnb$theta
qgamma((1:3)/4, shape = 1/v, scale = v)
```

    ## [1] 0.3116183 0.7149568 1.3850006

``` r
round(data.frame(p=coef(mp),q=coef(mq),nb=coef(mnb), se.p=se(mp),se.q=se(mq),se.nb=se(mnb)),4)
```

    ##                      p       q      nb   se.p   se.q  se.nb
    ## (Intercept)    12.6941 12.6941 15.9416 0.0424 4.0964 2.8558
    ## Gest_age       -0.1594 -0.1594 -0.2552 0.0013 0.1258 0.0850
    ## DIET_CODES_4F   1.5566  1.5566  1.4979 0.0110 1.0647 0.4799
    ## DIET_CODES_4FF  0.7888  0.7888  0.5061 0.0109 1.0548 0.5164

``` r
deviance(mnbg)
```

    ## [1] 52.59779

``` r
######################################
#Build the models for ARG abundance with overdispersion

# Gamma distrubuted GLM with log link
fit1_gamma<-glm(SUM2~DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit1_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.0877  -1.1826  -0.4111   0.3076   2.5556  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     -2.0636     0.5141  -4.014 0.000235 ***
    ## DIET_CODES_4F    1.2678     0.5721   2.216 0.032024 *  
    ## DIET_CODES_4FF   0.7040     0.5748   1.225 0.227339    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.321714)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 54.579  on 43  degrees of freedom
    ## AIC: -7.8464
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
summ(fit1_gamma, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 46
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(2) = 6.85, p = 0.08
    ## Pseudo-R² (Cragg-Uhler) = -0.58
    ## Pseudo-R² (McFadden) = -0.69
    ## AIC = -7.85, BIC = -0.53 
    ## 
    ## Standard errors: MLE
    ## ---------------------------------------------------------------
    ##                        exp(Est.)   2.5%   97.5%   t val.      p
    ## -------------------- ----------- ------ ------- -------- ------
    ## (Intercept)                 0.13   0.05    0.35    -4.01   0.00
    ## DIET_CODES_4F               3.55   1.16   10.90     2.22   0.03
    ## DIET_CODES_4FF              2.02   0.66    6.24     1.22   0.23
    ## ---------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.32

``` r
# Tukey
glht.mod <- glht(fit1_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## F - EBM == 0    1.2678     0.5721   2.216   0.0657 .
    ## FF - EBM == 0   0.7040     0.5748   1.225   0.4306  
    ## FF - F == 0    -0.5638     0.3592  -1.570   0.2520  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
fit2_gamma<-glm(SUM2~Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit2_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age + DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.9577  -0.9639  -0.1780   0.1088   2.0395  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)     7.44383    3.04181   2.447  0.01866 * 
    ## Gest_age       -0.29420    0.09058  -3.248  0.00229 **
    ## DIET_CODES_4F   1.59806    0.51110   3.127  0.00321 **
    ## DIET_CODES_4FF  0.25678    0.54994   0.467  0.64297   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.045734)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 46.436  on 42  degrees of freedom
    ## AIC: -14.521
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
summ(fit2_gamma, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 46
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(3) = 14.99, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = -1.24
    ## Pseudo-R² (McFadden) = -1.61
    ## AIC = -14.52, BIC = -5.38 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------------
    ##                        exp(Est.)   2.5%       97.5%   t val.      p
    ## -------------------- ----------- ------ ----------- -------- ------
    ## (Intercept)              1709.29   4.40   663750.82     2.45   0.02
    ## Gest_age                    0.75   0.62        0.89    -3.25   0.00
    ## DIET_CODES_4F               4.94   1.82       13.46     3.13   0.00
    ## DIET_CODES_4FF              1.29   0.44        3.80     0.47   0.64
    ## -------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.05

``` r
fit2_gamma_MGE<-glm(MGE_SUM~Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit2_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age + DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.9577  -0.9639  -0.1780   0.1088   2.0395  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)     7.44383    3.04181   2.447  0.01866 * 
    ## Gest_age       -0.29420    0.09058  -3.248  0.00229 **
    ## DIET_CODES_4F   1.59806    0.51110   3.127  0.00321 **
    ## DIET_CODES_4FF  0.25678    0.54994   0.467  0.64297   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.045734)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 46.436  on 42  degrees of freedom
    ## AIC: -14.521
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
summ(fit2_gamma, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 46
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(3) = 14.99, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = -1.24
    ## Pseudo-R² (McFadden) = -1.61
    ## AIC = -14.52, BIC = -5.38 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------------
    ##                        exp(Est.)   2.5%       97.5%   t val.      p
    ## -------------------- ----------- ------ ----------- -------- ------
    ## (Intercept)              1709.29   4.40   663750.82     2.45   0.02
    ## Gest_age                    0.75   0.62        0.89    -3.25   0.00
    ## DIET_CODES_4F               4.94   1.82       13.46     3.13   0.00
    ## DIET_CODES_4FF              1.29   0.44        3.80     0.47   0.64
    ## -------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.05

``` r
glht.mod <- glht(fit2_gamma_MGE, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## F - EBM == 0   1.18476    0.56939   2.081   0.0921 .
    ## FF - EBM == 0  0.08353    0.61266   0.136   0.9897  
    ## FF - F == 0   -1.10123    0.45215  -2.436   0.0385 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
qchisq(0.95, df.residual(fit2_gamma))
```

    ## [1] 58.12404

``` r
deviance(fit2_gamma)
```

    ## [1] 46.43596

``` r
cor(model.matrix(fit2_gamma)[,-1])
```

    ##                  Gest_age DIET_CODES_4F DIET_CODES_4FF
    ## Gest_age        1.0000000      0.541648     -0.6195105
    ## DIET_CODES_4F   0.5416480      1.000000     -0.8038370
    ## DIET_CODES_4FF -0.6195105     -0.803837      1.0000000

``` r
glht.mod <- glht(fit2_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##               Estimate Std. Error z value Pr(>|z|)   
    ## F - EBM == 0    1.5981     0.5111   3.127  0.00491 **
    ## FF - EBM == 0   0.2568     0.5499   0.467  0.88542   
    ## FF - F == 0    -1.3413     0.4059  -3.305  0.00273 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
fit2_gamma_ntw<-glm(SUM2~Gest_age+DIET_CODES_4, data=df[df$Twin=="n",], family=Gamma(link="log"))
summary(fit2_gamma_ntw)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age + DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df[df$Twin == "n", ])
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.83867  -0.68360  -0.17033   0.07316   1.83688  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)     6.99313    2.95380   2.368  0.02325 * 
    ## Gest_age       -0.28067    0.08799  -3.190  0.00290 **
    ## DIET_CODES_4F   1.59751    0.49289   3.241  0.00252 **
    ## DIET_CODES_4FF  0.42078    0.53769   0.783  0.43886   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.9528069)
    ## 
    ##     Null deviance: 47.517  on 40  degrees of freedom
    ## Residual deviance: 35.133  on 37  degrees of freedom
    ## AIC: -8.3798
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
summ(fit2_gamma_ntw)
```

    ## MODEL INFO:
    ## Observations: 41
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(3) = 12.38, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = -2.81
    ## Pseudo-R² (McFadden) = -3.49
    ## AIC = -8.38, BIC = 0.19 
    ## 
    ## Standard errors: MLE
    ## ---------------------------------------------------
    ##                         Est.   S.E.   t val.      p
    ## -------------------- ------- ------ -------- ------
    ## (Intercept)             6.99   2.95     2.37   0.02
    ## Gest_age               -0.28   0.09    -3.19   0.00
    ## DIET_CODES_4F           1.60   0.49     3.24   0.00
    ## DIET_CODES_4FF          0.42   0.54     0.78   0.44
    ## ---------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 0.95

``` r
fit2_gamma_nof<-glm(SUM2~Gest_age+DIET_CODES_4, data=df[df$Formula=="n",], family=Gamma(link="log"))
summary(fit2_gamma_nof)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age + DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df[df$Formula == "n", ])
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.88811  -0.98299  -0.06715   0.05294   2.17681  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)      5.4181     3.4182   1.585   0.1272  
    ## Gest_age        -0.2331     0.1019  -2.288   0.0321 *
    ## DIET_CODES_4FF   0.3813     0.6000   0.636   0.5316  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.232576)
    ## 
    ##     Null deviance: 29.467  on 24  degrees of freedom
    ## Residual deviance: 22.852  on 22  degrees of freedom
    ## AIC: -22.354
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
summ(fit2_gamma_nof)
```

    ## MODEL INFO:
    ## Observations: 25
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(2) = 6.61, p = 0.07
    ## Pseudo-R² (Cragg-Uhler) = -0.17
    ## Pseudo-R² (McFadden) = -0.32
    ## AIC = -22.35, BIC = -17.48 
    ## 
    ## Standard errors: MLE
    ## ---------------------------------------------------
    ##                         Est.   S.E.   t val.      p
    ## -------------------- ------- ------ -------- ------
    ## (Intercept)             5.42   3.42     1.59   0.13
    ## Gest_age               -0.23   0.10    -2.29   0.03
    ## DIET_CODES_4FF          0.38   0.60     0.64   0.53
    ## ---------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.23

``` r
library(pscl)
```

    ## Classes and Methods for R developed in the
    ## Political Science Computational Laboratory
    ## Department of Political Science
    ## Stanford University
    ## Simon Jackman
    ## hurdle and zeroinfl functions by Achim Zeileis

``` r
pchisq(summary(fit2_gamma)$deviance, 
           summary(fit2_gamma)$df.residual
           )
```

    ## [1] 0.7054517

``` r
fit3_gamma<-glm(SUM2~Gest_age*DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit3_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age * DIET_CODES_4, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8934  -1.0259  -0.1318   0.2614   2.2280  
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)               6.6747     4.8426   1.378  0.17576   
    ## Gest_age                 -0.2711     0.1452  -1.866  0.06932 . 
    ## DIET_CODES_4F            36.4602    12.3975   2.941  0.00542 **
    ## DIET_CODES_4FF           -1.6346     6.0643  -0.270  0.78889   
    ## Gest_age:DIET_CODES_4F   -1.0378     0.3682  -2.819  0.00746 **
    ## Gest_age:DIET_CODES_4FF   0.0625     0.1869   0.334  0.73987   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.9929224)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 40.585  on 40  degrees of freedom
    ## AIC: -17.625
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
summ(fit3_gamma)
```

    ## MODEL INFO:
    ## Observations: 46
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(5) = 20.84, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = -1.69
    ## Pseudo-R² (McFadden) = -2.37
    ## AIC = -17.63, BIC = -4.82 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------
    ##                                  Est.    S.E.   t val.      p
    ## ----------------------------- ------- ------- -------- ------
    ## (Intercept)                      6.67    4.84     1.38   0.18
    ## Gest_age                        -0.27    0.15    -1.87   0.07
    ## DIET_CODES_4F                   36.46   12.40     2.94   0.01
    ## DIET_CODES_4FF                  -1.63    6.06    -0.27   0.79
    ## Gest_age:DIET_CODES_4F          -1.04    0.37    -2.82   0.01
    ## Gest_age:DIET_CODES_4FF          0.06    0.19     0.33   0.74
    ## -------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 0.99

``` r
pR2(fit3_gamma)
```

    ##       llh   llhNull        G2  McFadden      r2ML      r2CU 
    ## 15.812660  4.694419 22.236482 -2.368396  0.383318 -1.692898

``` r
temp<-cbind(df, 
      Mean = predict(fit3_gamma, newdata=df, type="response"), 
      SE = predict(fit3_gamma, newdata=df, type="response", se.fit=T)$se.fit
      )

fit4<-glm.nb((SUM)~Inf_AB+Gest_age*DIET_CODES_4+offset(log(SSU_counts)), data=df, link = log)
summary(fit4)
```

    ## 
    ## Call:
    ## glm.nb(formula = (SUM) ~ Inf_AB + Gest_age * DIET_CODES_4 + offset(log(SSU_counts)), 
    ##     data = df, link = log, init.theta = 1.423143806)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.0821  -1.0690  -0.3445   0.3504   2.8290  
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              4.50194    4.07741   1.104 0.269543    
    ## Inf_ABY                 -0.24464    0.27352  -0.894 0.371089    
    ## Gest_age                -0.22530    0.12236  -1.841 0.065589 .  
    ## DIET_CODES_4F           37.41519   10.55918   3.543 0.000395 ***
    ## DIET_CODES_4FF           0.30761    5.17002   0.059 0.952554    
    ## Gest_age:DIET_CODES_4F  -1.06085    0.31387  -3.380 0.000725 ***
    ## Gest_age:DIET_CODES_4FF  0.01232    0.15965   0.077 0.938491    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(1.4231) family taken to be 1)
    ## 
    ##     Null deviance: 79.335  on 45  degrees of freedom
    ## Residual deviance: 51.150  on 39  degrees of freedom
    ## AIC: 879.92
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  1.423 
    ##           Std. Err.:  0.269 
    ## 
    ##  2 x log-likelihood:  -863.925

``` r
fit4_gamma<-glm(SUM2~Age_sampling+Inf_AB+Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit4_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Age_sampling + Inf_AB + Gest_age + DIET_CODES_4, 
    ##     family = Gamma(link = "log"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.9304  -0.9501  -0.2020   0.1301   1.8426  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)    10.20118    4.17489   2.443  0.01905 * 
    ## Age_sampling   -0.05226    0.05331  -0.980  0.33285   
    ## Inf_ABY         0.15819    0.38317   0.413  0.68193   
    ## Gest_age       -0.35263    0.10805  -3.264  0.00226 **
    ## DIET_CODES_4F   1.49314    0.53657   2.783  0.00819 **
    ## DIET_CODES_4FF  0.08835    0.59694   0.148  0.88308   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.041566)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 45.629  on 40  degrees of freedom
    ## AIC: -11.452
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
summ(fit4_gamma, exp = TRUE)
```

    ## MODEL INFO:
    ## Observations: 46
    ## Dependent Variable: SUM2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(5) = 15.80, p = 0.01
    ## Pseudo-R² (Cragg-Uhler) = -1.30
    ## Pseudo-R² (McFadden) = -1.71
    ## AIC = -11.45, BIC = 1.35 
    ## 
    ## Standard errors: MLE
    ## ---------------------------------------------------------------------
    ##                        exp(Est.)   2.5%         97.5%   t val.      p
    ## -------------------- ----------- ------ ------------- -------- ------
    ## (Intercept)             26934.97   7.53   96379685.59     2.44   0.02
    ## Age_sampling                0.95   0.85          1.05    -0.98   0.33
    ## Inf_ABY                     1.17   0.55          2.48     0.41   0.68
    ## Gest_age                    0.70   0.57          0.87    -3.26   0.00
    ## DIET_CODES_4F               4.45   1.56         12.74     2.78   0.01
    ## DIET_CODES_4FF              1.09   0.34          3.52     0.15   0.88
    ## ---------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.04

``` r
fit5_gamma <- glm(SUM2~Age_sampling+Inf_AB+Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit5_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Age_sampling + Inf_AB + Gest_age + DIET_CODES_4, 
    ##     family = Gamma(link = "log"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.9304  -0.9501  -0.2020   0.1301   1.8426  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)    10.20118    4.17489   2.443  0.01905 * 
    ## Age_sampling   -0.05226    0.05331  -0.980  0.33285   
    ## Inf_ABY         0.15819    0.38317   0.413  0.68193   
    ## Gest_age       -0.35263    0.10805  -3.264  0.00226 **
    ## DIET_CODES_4F   1.49314    0.53657   2.783  0.00819 **
    ## DIET_CODES_4FF  0.08835    0.59694   0.148  0.88308   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.041566)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 45.629  on 40  degrees of freedom
    ## AIC: -11.452
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
fit6_gamma <- glm(SUM2~Age_sampling+Mat_AB+Gest_age+DIET_CODES_4, data=df, family=Gamma(link="log"))
summary(fit6_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Age_sampling + Mat_AB + Gest_age + DIET_CODES_4, 
    ##     family = Gamma(link = "log"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.0258  -0.8093  -0.2200   0.2416   1.6592  
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    10.98835    3.81313   2.882 0.006332 ** 
    ## Age_sampling   -0.04125    0.04186  -0.985 0.330351    
    ## Mat_ABY         0.50262    0.30245   1.662 0.104365    
    ## Gest_age       -0.38341    0.09998  -3.835 0.000436 ***
    ## DIET_CODES_4F   1.41611    0.49818   2.843 0.007013 ** 
    ## DIET_CODES_4FF -0.20030    0.57757  -0.347 0.730557    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.8844819)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 43.462  on 40  degrees of freedom
    ## AIC: -14.026
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
# Perform chisquared test

anova_NEC<-anova(fit1_gamma, fit2_gamma, fit3_gamma, fit4_gamma, test="Chisq")


# Use Tukey's post hoc test for the model 2
glht.mod <- glht(fit2_gamma, mcp(DIET_CODES_4="Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##               Estimate Std. Error z value Pr(>|z|)   
    ## F - EBM == 0    1.5981     0.5111   3.127  0.00486 **
    ## FF - EBM == 0   0.2568     0.5499   0.467  0.88542   
    ## FF - F == 0    -1.3413     0.4059  -3.305  0.00268 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
fit5_gamma<-glm(SUM2~Gest_age+DIET_CODES_1, data=df, family=Gamma(link="log"))
summary(fit5_gamma)
```

    ## 
    ## Call:
    ## glm(formula = SUM2 ~ Gest_age + DIET_CODES_1, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.7867  -1.0205  -0.1185   0.3444   1.2892  
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                           8.05097    2.97736   2.704  0.01062
    ## Gest_age                             -0.30048    0.09344  -3.216  0.00285
    ## DIET_CODES_1exclusively breastfed    -0.39758    1.00112  -0.397  0.69375
    ## DIET_CODES_1formula                   0.63257    1.29230   0.489  0.62764
    ## DIET_CODES_1mom +HMF                 -0.52349    0.95657  -0.547  0.58777
    ## DIET_CODES_1mom +Prolacta             0.30236    0.98798   0.306  0.76144
    ## DIET_CODES_1predom donor+mom+HMF      0.06432    1.00713   0.064  0.94945
    ## DIET_CODES_1predom formula +mom       1.13658    0.97161   1.170  0.25022
    ## DIET_CODES_1predom formula+mom+donor -1.23860    1.12774  -1.098  0.27979
    ## DIET_CODES_1predom mom +donor+HMF    -0.75138    1.02003  -0.737  0.46640
    ## DIET_CODES_1predom mom +formula       1.60756    0.99510   1.615  0.11545
    ## DIET_CODES_1predom mom+donor+HMF     -0.29460    1.12034  -0.263  0.79417
    ##                                        
    ## (Intercept)                          * 
    ## Gest_age                             **
    ## DIET_CODES_1exclusively breastfed      
    ## DIET_CODES_1formula                    
    ## DIET_CODES_1mom +HMF                   
    ## DIET_CODES_1mom +Prolacta              
    ## DIET_CODES_1predom donor+mom+HMF       
    ## DIET_CODES_1predom formula +mom        
    ## DIET_CODES_1predom formula+mom+donor   
    ## DIET_CODES_1predom mom +donor+HMF      
    ## DIET_CODES_1predom mom +formula        
    ## DIET_CODES_1predom mom+donor+HMF       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.7793765)
    ## 
    ##     Null deviance: 61.424  on 45  degrees of freedom
    ## Residual deviance: 36.017  on 34  degrees of freedom
    ## AIC: -11.839
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
glht.mod <- glht(fit5_gamma, mcp(DIET_CODES_1="Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##                                                       Estimate Std. Error
    ## exclusively breastfed - donor +HMF == 0               -0.39758    1.00112
    ## formula - donor +HMF == 0                              0.63257    1.29230
    ## mom +HMF - donor +HMF == 0                            -0.52349    0.95657
    ## mom +Prolacta - donor +HMF == 0                        0.30236    0.98798
    ## predom donor+mom+HMF - donor +HMF == 0                 0.06432    1.00713
    ## predom formula +mom - donor +HMF == 0                  1.13658    0.97161
    ## predom formula+mom+donor - donor +HMF == 0            -1.23860    1.12774
    ## predom mom +donor+HMF - donor +HMF == 0               -0.75138    1.02003
    ## predom mom +formula - donor +HMF == 0                  1.60756    0.99510
    ## predom mom+donor+HMF - donor +HMF == 0                -0.29460    1.12034
    ## formula - exclusively breastfed == 0                   1.03015    0.96997
    ## mom +HMF - exclusively breastfed == 0                 -0.12592    0.63070
    ## mom +Prolacta - exclusively breastfed == 0             0.69994    0.63017
    ## predom donor+mom+HMF - exclusively breastfed == 0      0.46190    0.59511
    ## predom formula +mom - exclusively breastfed == 0       1.53416    0.47335
    ## predom formula+mom+donor - exclusively breastfed == 0 -0.84102    0.74119
    ## predom mom +donor+HMF - exclusively breastfed == 0    -0.35381    0.70886
    ## predom mom +formula - exclusively breastfed == 0       2.00514    0.53519
    ## predom mom+donor+HMF - exclusively breastfed == 0      0.10297    0.73943
    ## mom +HMF - formula == 0                               -1.15606    1.03774
    ## mom +Prolacta - formula == 0                          -0.33021    1.02879
    ## predom donor+mom+HMF - formula == 0                   -0.56825    0.99600
    ## predom formula +mom - formula == 0                     0.50402    0.91904
    ## predom formula+mom+donor - formula == 0               -1.87117    1.08131
    ## predom mom +donor+HMF - formula == 0                  -1.38395    1.08427
    ## predom mom +formula - formula == 0                     0.97500    0.95482
    ## predom mom+donor+HMF - formula == 0                   -0.92717    1.08198
    ## mom +Prolacta - mom +HMF == 0                          0.82585    0.58221
    ## predom donor+mom+HMF - mom +HMF == 0                   0.58781    0.63320
    ## predom formula +mom - mom +HMF == 0                    1.66008    0.59008
    ## predom formula+mom+donor - mom +HMF == 0              -0.71510    0.82261
    ## predom mom +donor+HMF - mom +HMF == 0                 -0.22789    0.62553
    ## predom mom +formula - mom +HMF == 0                    2.13106    0.62423
    ## predom mom+donor+HMF - mom +HMF == 0                   0.22889    0.80990
    ## predom donor+mom+HMF - mom +Prolacta == 0             -0.23804    0.64363
    ## predom formula +mom - mom +Prolacta == 0               0.83423    0.57788
    ## predom formula+mom+donor - mom +Prolacta == 0         -1.54096    0.81320
    ## predom mom +donor+HMF - mom +Prolacta == 0            -1.05374    0.67891
    ## predom mom +formula - mom +Prolacta == 0               1.30520    0.61876
    ## predom mom+donor+HMF - mom +Prolacta == 0             -0.59696    0.80437
    ## predom formula +mom - predom donor+mom+HMF == 0        1.07227    0.52263
    ## predom formula+mom+donor - predom donor+mom+HMF == 0  -1.30292    0.77395
    ## predom mom +donor+HMF - predom donor+mom+HMF == 0     -0.81570    0.71438
    ## predom mom +formula - predom donor+mom+HMF == 0        1.54325    0.57607
    ## predom mom+donor+HMF - predom donor+mom+HMF == 0      -0.35892    0.77021
    ## predom formula+mom+donor - predom formula +mom == 0   -2.37518    0.67428
    ## predom mom +donor+HMF - predom formula +mom == 0      -1.88797    0.66959
    ## predom mom +formula - predom formula +mom == 0         0.47098    0.44251
    ## predom mom+donor+HMF - predom formula +mom == 0       -1.43119    0.67464
    ## predom mom +donor+HMF - predom formula+mom+donor == 0  0.48722    0.88116
    ## predom mom +formula - predom formula+mom+donor == 0    2.84616    0.72172
    ## predom mom+donor+HMF - predom formula+mom+donor == 0   0.94400    0.88324
    ## predom mom +formula - predom mom +donor+HMF == 0       2.35895    0.70165
    ## predom mom+donor+HMF - predom mom +donor+HMF == 0      0.45678    0.87056
    ## predom mom+donor+HMF - predom mom +formula == 0       -1.90217    0.72088
    ##                                                       z value Pr(>|z|)   
    ## exclusively breastfed - donor +HMF == 0                -0.397   1.0000   
    ## formula - donor +HMF == 0                               0.489   1.0000   
    ## mom +HMF - donor +HMF == 0                             -0.547   1.0000   
    ## mom +Prolacta - donor +HMF == 0                         0.306   1.0000   
    ## predom donor+mom+HMF - donor +HMF == 0                  0.064   1.0000   
    ## predom formula +mom - donor +HMF == 0                   1.170   0.9829   
    ## predom formula+mom+donor - donor +HMF == 0             -1.098   0.9894   
    ## predom mom +donor+HMF - donor +HMF == 0                -0.737   0.9996   
    ## predom mom +formula - donor +HMF == 0                   1.615   0.8599   
    ## predom mom+donor+HMF - donor +HMF == 0                 -0.263   1.0000   
    ## formula - exclusively breastfed == 0                    1.062   0.9918   
    ## mom +HMF - exclusively breastfed == 0                  -0.200   1.0000   
    ## mom +Prolacta - exclusively breastfed == 0              1.111   0.9884   
    ## predom donor+mom+HMF - exclusively breastfed == 0       0.776   0.9994   
    ## predom formula +mom - exclusively breastfed == 0        3.241   0.0409 * 
    ## predom formula+mom+donor - exclusively breastfed == 0  -1.135   0.9864   
    ## predom mom +donor+HMF - exclusively breastfed == 0     -0.499   1.0000   
    ## predom mom +formula - exclusively breastfed == 0        3.747    <0.01 **
    ## predom mom+donor+HMF - exclusively breastfed == 0       0.139   1.0000   
    ## mom +HMF - formula == 0                                -1.114   0.9881   
    ## mom +Prolacta - formula == 0                           -0.321   1.0000   
    ## predom donor+mom+HMF - formula == 0                    -0.571   1.0000   
    ## predom formula +mom - formula == 0                      0.548   1.0000   
    ## predom formula+mom+donor - formula == 0                -1.730   0.7991   
    ## predom mom +donor+HMF - formula == 0                   -1.276   0.9679   
    ## predom mom +formula - formula == 0                      1.021   0.9940   
    ## predom mom+donor+HMF - formula == 0                    -0.857   0.9986   
    ## mom +Prolacta - mom +HMF == 0                           1.418   0.9352   
    ## predom donor+mom+HMF - mom +HMF == 0                    0.928   0.9972   
    ## predom formula +mom - mom +HMF == 0                     2.813   0.1369   
    ## predom formula+mom+donor - mom +HMF == 0               -0.869   0.9984   
    ## predom mom +donor+HMF - mom +HMF == 0                  -0.364   1.0000   
    ## predom mom +formula - mom +HMF == 0                     3.414   0.0234 * 
    ## predom mom+donor+HMF - mom +HMF == 0                    0.283   1.0000   
    ## predom donor+mom+HMF - mom +Prolacta == 0              -0.370   1.0000   
    ## predom formula +mom - mom +Prolacta == 0                1.444   0.9277   
    ## predom formula+mom+donor - mom +Prolacta == 0          -1.895   0.6948   
    ## predom mom +donor+HMF - mom +Prolacta == 0             -1.552   0.8882   
    ## predom mom +formula - mom +Prolacta == 0                2.109   0.5406   
    ## predom mom+donor+HMF - mom +Prolacta == 0              -0.742   0.9996   
    ## predom formula +mom - predom donor+mom+HMF == 0         2.052   0.5824   
    ## predom formula+mom+donor - predom donor+mom+HMF == 0   -1.683   0.8255   
    ## predom mom +donor+HMF - predom donor+mom+HMF == 0      -1.142   0.9857   
    ## predom mom +formula - predom donor+mom+HMF == 0         2.679   0.1892   
    ## predom mom+donor+HMF - predom donor+mom+HMF == 0       -0.466   1.0000   
    ## predom formula+mom+donor - predom formula +mom == 0    -3.523   0.0165 * 
    ## predom mom +donor+HMF - predom formula +mom == 0       -2.820   0.1351   
    ## predom mom +formula - predom formula +mom == 0          1.064   0.9916   
    ## predom mom+donor+HMF - predom formula +mom == 0        -2.121   0.5310   
    ## predom mom +donor+HMF - predom formula+mom+donor == 0   0.553   1.0000   
    ## predom mom +formula - predom formula+mom+donor == 0     3.944    <0.01 **
    ## predom mom+donor+HMF - predom formula+mom+donor == 0    1.069   0.9914   
    ## predom mom +formula - predom mom +donor+HMF == 0        3.362   0.0280 * 
    ## predom mom+donor+HMF - predom mom +donor+HMF == 0       0.525   1.0000   
    ## predom mom+donor+HMF - predom mom +formula == 0        -2.639   0.2088   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
library(ggsignif)

annotation_df <- data_frame(start=c("FF", "EBM", "EBM"), 
                            end=c("F", "F", "FF"),
                            y=c(0.5, 0.75, 1),
                            label=c("**", "**", "ns"))
```

    ## Warning: `data_frame()` is deprecated, use `tibble()`.
    ## This warning is displayed once per session.

``` r
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
```

    ## Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position

``` r
ARG.sum.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Plot efect of gestational age

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
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
# Build  model for MGE abundace

fit1<-glm((MGE_SUM)/mean16S_counts~Formula, data=df, family="Gamma"(link="log"))
summary(fit1)
```

    ## 
    ## Call:
    ## glm(formula = (MGE_SUM)/mean16S_counts ~ Formula, family = Gamma(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.4198  -1.4910  -0.5309   0.3538   1.8381  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.8390     0.2259   3.714 0.000573 ***
    ## Formulay      0.7292     0.3344   2.181 0.034593 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.275998)
    ## 
    ##     Null deviance: 76.309  on 45  degrees of freedom
    ## Residual deviance: 70.247  on 44  degrees of freedom
    ## AIC: 204.66
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
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
```

    ## Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position

``` r
MGE.sum.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
# Build  model for Enterobacteriaceae abundace
fit1<-glm(ENTERO/100~Gest_age+Formula, data=df, family="quasibinomial")
summary(fit1)
```

    ## 
    ## Call:
    ## glm(formula = ENTERO/100 ~ Gest_age + Formula, family = "quasibinomial", 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.4263  -0.6941   0.3757   0.6658   0.8891  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   8.0253     4.2044   1.909   0.0630 .
    ## Gest_age     -0.2543     0.1335  -1.905   0.0634 .
    ## Formulay      1.1397     0.5462   2.087   0.0429 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 0.5266167)
    ## 
    ##     Null deviance: 29.760  on 45  degrees of freedom
    ## Residual deviance: 26.803  on 43  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 3

``` r
# Build model for Gammaproteobacteria abundance
fit1_gamma<-glm(GAMMA/100~Gest_age+Formula, data=df, family="quasibinomial")
summary(fit1_gamma)
```

    ## 
    ## Call:
    ## glm(formula = GAMMA/100 ~ Gest_age + Formula, family = "quasibinomial", 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.4225  -0.7253   0.3528   0.6456   1.0805  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   8.2846     4.2231   1.962   0.0563 .
    ## Gest_age     -0.2579     0.1337  -1.929   0.0604 .
    ## Formulay      1.0437     0.5410   1.929   0.0603 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 0.5226474)
    ## 
    ##     Null deviance: 28.693  on 45  degrees of freedom
    ## Residual deviance: 25.959  on 43  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 3

``` r
# Make cowplot
cowplot:::plot_grid(ARG.sum.plot, MGE.sum.plot, labels="auto", align="h", nrow=1)
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

# DeSEQ

``` r
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

![](Supplementary_software_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Mantel’s test for correlations between taxa, ARGs and MGEs

``` r
#Do a mantel test for ARG and MGE distance matrixes

MGE_dist<-vegdist(t(as.matrix(otu_table(MGE_PHY))), method="horn")
ARG_dist<-vegdist(t(as.matrix(otu_table(ARG_PHY_res2))), method="horn")
mantel(ARG_dist, MGE_dist, method="kendall")
```

    ## 
    ## Mantel statistic based on Kendall's rank correlation tau 
    ## 
    ## Call:
    ## mantel(xdis = ARG_dist, ydis = MGE_dist, method = "kendall") 
    ## 
    ## Mantel statistic r:  0.23 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0675 0.0817 0.0990 0.1135 
    ## Permutation: free
    ## Number of permutations: 999

``` r
PHY_dist<-vegdist(t(as.matrix(otu_table(PHY_SP))), method="horn")
ARG_dist<-vegdist(t(as.matrix(otu_table(ARG_PHY_res2))), method="horn")
mantel(ARG_dist, PHY_dist, method="kendall")
```

    ## 
    ## Mantel statistic based on Kendall's rank correlation tau 
    ## 
    ## Call:
    ## mantel(xdis = ARG_dist, ydis = PHY_dist, method = "kendall") 
    ## 
    ## Mantel statistic r: 0.208 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0436 0.0568 0.0686 0.0872 
    ## Permutation: free
    ## Number of permutations: 999

``` r
PHY_dist<-vegdist(t(as.matrix(otu_table(PHY_SP))), method="horn")
MGE_dist<-vegdist(t(as.matrix(otu_table(MGE_PHY))), method="horn")
mantel(MGE_dist, PHY_dist, method="kendall")
```

    ## 
    ## Mantel statistic based on Kendall's rank correlation tau 
    ## 
    ## Call:
    ## mantel(xdis = MGE_dist, ydis = PHY_dist, method = "kendall") 
    ## 
    ## Mantel statistic r: 0.2635 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0504 0.0608 0.0764 0.0955 
    ## Permutation: free
    ## Number of permutations: 999

# Similarity between twins

``` r
# Similarity between twins using vegdist

# Read in function for testing significances
lmp = function(y,x,n.perms = 9999){
   test = numeric(n.perms)
   test[1] = abs(summary(lm(y~x))$coefficients[2,3])
   for(ii in 2:n.perms){
       test[ii] = abs(summary(lm(sample(y)~x))$coefficients[2,3])
   }
   p.value = sum(test>=test[1])/n.perms
   return(p.value)
}

# Species
# make a distance matrix for samples
otu <- sqrt(t(otu_table(PHY_SP)))
y <- vegdist(otu, upper=TRUE, diag=TRUE, method="horn")
distance <- as.vector(y)

# make a matrix for family aka PAIR
x = matrix(NA,46,46)
Fvec <- sample_data(PHY_SP)$Pair
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Fvec[n] == Fvec[m])
    x[n,m] = f
  }
}
x <- t(x)
family_vec <- x[lower.tri(x,diag=F)]

# and the same for DIET
x = matrix(NA,46,46)
Yvec <- sample_data(PHY_SP)$Formula
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Yvec[n] == Yvec[m])
    x[n,m] = f
  }
}
x <- t(x)
type_vec <- x[lower.tri(x,diag=F)]


# First is for family and second is for type. 
zz <-data.frame(distance, family_vec, type_vec)
zz$comb <- paste(zz$family_vec, zz$type_vec, sep="-")
zz<-zz[grep("TRUE-FALSE", zz$comb, invert = TRUE),]


a1<-ggplot(zz, aes(x=distance, fill=factor(comb))) + geom_density(alpha=0.5) + theme(legend.justification=c(0.05,1), legend.position=c(0.05,1), legend.text= element_text(size=rel(0.8))) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"),  "Species")  +ggtitle("Species") + ylim(0,20)


# Test significance

lmp(zz$distance, zz$family_vec)
```

    ## [1] 0.00410041

``` r
#######
# Remove the FALSE-TRUE, or different family same type, this is significant 0.0001 (Difference between twins with same diet and non twins dif diet)
zz_temp<-zz[grep("FALSE-TRUE", zz$comb, invert = TRUE),] 

# Test significance
lmp(zz_temp$distance, zz_temp$comb)
```

    ## [1] 0.00420042

``` r
# Remove the FALSE-FALSE,
zz_temp<-zz[grep("FALSE-FALSE", zz$comb, invert = TRUE),] 

# Test significance, p=0.0002, there is difference. Twins more similar to each other than to unrelated infants with same diet
lmp(zz_temp$distance, zz_temp$comb)
```

    ## [1] 0.00390039

``` r
# Remove TRUE-FALSE
zz_temp<-zz[grep("TRUE-TRUE", zz$comb, invert = TRUE),] 

# Test significance, no difference between unrelated infants with different diets.
lmp(zz_temp$distance, zz_temp$comb)
```

    ## [1] 0.7381738

``` r
###########

# ARGs
# make a distance matrix for samples
otu <- (t(otu_table(ARG_PHY_res2)))
y <- vegdist(otu, upper=TRUE, diag=TRUE, method="horn")
distance <- as.vector(y)

# Make a matrix for family aka PAIR
x = matrix(NA,46,46)
Fvec <- sample_data(ARG_PHY_res2)$Pair
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Fvec[n] == Fvec[m])
    x[n,m] = f
  }
}
x <- t(x)
family_vec <- x[lower.tri(x,diag=F)]

# and the same for diet
x = matrix(NA,46,46)
Yvec <- sample_data(ARG_PHY_res2)$Formula
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Yvec[n] == Yvec[m])
    x[n,m] = f
  }
}
x <- t(x)
type_vec <- x[lower.tri(x,diag=F)]


# First is for family and second is for type. 
zz <-data.frame(distance, family_vec, type_vec)
zz$comb <- paste(zz$family_vec, zz$type_vec, sep="-")
zz<-zz[grep("TRUE-FALSE", zz$comb, invert = TRUE),]


b1<-ggplot(zz, aes(x=distance, fill=factor(comb))) + geom_density(alpha=0.5) + theme(legend.justification=c(0.05,1), legend.position=c(0.05,1), legend.text= element_text(size=rel(0.8))) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"),  "ARGs") + xlab("Between-sample dissimilarity") +ggtitle("ARGs") + ylim(0,20)

# Test significance
lmp(zz$distance, zz$family_vec)
```

    ## [1] 0.00150015

``` r
# Remove TRUE-FALSE, or same family different type, this is significant p=0.0001, dif. fam dif type vs. dif fam same type, type is significant
zz_temp<-zz[grep("TRUE-TRUE", zz$comb, invert = TRUE),] 

 # Test significance, no difference between unrelated infants with different diets.
lmp(zz_temp$distance, zz_temp$comb)
```

    ## [1] 0.9407941

``` r
# MGEs
# make a distance matrix for samples
otu <- (t(otu_table(MGE_PHY)))
y <- vegdist(otu, upper=TRUE, diag=TRUE, method="horn")
distance <- as.vector(y)

# make a matrix for family aka PAIR
x = matrix(NA,46,46)
Fvec <- sample_data(MGE_PHY)$Pair
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Fvec[n] == Fvec[m])
    x[n,m] = f
  }
}
x <- t(x)
family_vec <- x[lower.tri(x,diag=F)]


# and the same for diet
x = matrix(NA,46,46)
Yvec <- sample_data(MGE_PHY)$Formula
for(n in 1:(ncol(x)-1)) {
  for(m in (n+1):ncol(x)){
    f <- eval(Yvec[n] == Yvec[m])
    x[n,m] = f
  }
}
x <- t(x)
type_vec <- x[lower.tri(x,diag=F)]

# First is for family and second is for type. 
zz <-data.frame(distance, family_vec, type_vec)
zz$comb <- paste(zz$family_vec, zz$type_vec, sep="-")
zz<-zz[grep("TRUE-FALSE", zz$comb, invert = TRUE),]

# Density plot
c1<-ggplot(zz, aes(x=distance, fill=factor(comb))) + geom_density(alpha=0.5) + theme(legend.justification=c(0.05,1), legend.position=c(0.05,1), legend.text= element_text(size=rel(0.8))) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"),  "MGEs") + xlab("Between-sample dissimilarity") +ggtitle("MGEs") + ylim(0,20)


# Test significance
lmp(zz$distance, zz$family_vec)
```

    ## [1] 0.00010001

``` r
cowplot::plot_grid(a1, b1, c1, labels="auto")
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

# Meta-analysis

Phyloseq objects of the meta-analysi dataset were created in similar way
as for the 46-infant cohort.

``` r
# Load meta-analysis Phyloseq object with ARG mapping results

merged_phyloseq <- readRDS("merged_phyloseq.rds")

#Make data frame

df_meta<-data.frame(SUM=sample_sums(merged_phyloseq),
                    Formula=as.factor(sample_data(merged_phyloseq)$Formula2), Antibiotics=sample_data(merged_phyloseq)$Antibiotics, Delivery_mode=sample_data(merged_phyloseq)$Delivery_mode,
                    Fortifier=sample_data(merged_phyloseq)$Fortifier,
                    Age=sample_data(merged_phyloseq)$Age,
                    Gest_age=sample_data(merged_phyloseq)$Gest_age,
                    Study=sample_data(merged_phyloseq)$Study,
                    Sample=sample_names(merged_phyloseq),
                    Donor=sample_data(merged_phyloseq)$Donor_Milk,
                    Country=sample_data(merged_phyloseq)$Country,
                    SSU_counts=sample_data(merged_phyloseq)$SSU_counts,
                    Twin=sample_data(merged_phyloseq)$Twin)

#Save df for fortifier analysis
df_meta_fort<-df_meta[df_meta$Study%in%c("NEC", "Gibson"),]

df_meta_fort<-df_meta_fort[df_meta_fort$Formula=="FALSE",]

#Exclude fortifier samples that don't have formula as well
df_meta<-df_meta[!(df_meta$Fortifier=="TRUE"&df_meta$Formula=="FALSE"),]
```

\#GLMs for meta

``` r
# Build models
df_meta$Formula<-as.logical(df_meta$Formula)
df_meta$Fortifier<-as.logical(df_meta$Fortifier)

#For analysis excluding twins
# df_meta<-df_meta[!(df_meta$Twin=="y"),]

fit_meta1<-glm(SUM~Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta1, exp=TRUE, digits=4)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(1) = 13.1729, p = 0.0109
    ## Pseudo-R² (Cragg-Uhler) = 0.1475
    ## Pseudo-R² (McFadden) = 0.1269
    ## AIC = 70.1275, BIC = 80.1111 
    ## 
    ## Standard errors: MLE
    ## ------------------------------------------------------------------
    ##                     exp(Est.)     2.5%    97.5%    t val.        p
    ## ----------------- ----------- -------- -------- --------- --------
    ## (Intercept)            0.3382   0.2539   0.4506   -7.4095   0.0000
    ## FormulaTRUE            1.6709   1.1306   2.4694    2.5759   0.0107
    ## ------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 2.0332

``` r
summary(fit_meta1)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ Formula, family = Gamma(link = "log"), data = df_meta)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.5406  -1.1837  -0.4736   0.0627   3.5449  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -1.0840     0.1463  -7.410  3.3e-12 ***
    ## FormulaTRUE   0.5134     0.1993   2.576   0.0107 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 2.033232)
    ## 
    ##     Null deviance: 369.30  on 205  degrees of freedom
    ## Residual deviance: 356.12  on 204  degrees of freedom
    ## AIC: 70.128
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
fit_meta2<-glm(SUM~Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta2, exp=TRUE, digits=4)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(2) = 34.4839, p = 0.0000
    ## Pseudo-R² (Cragg-Uhler) = 0.3817
    ## Pseudo-R² (McFadden) = 0.3410
    ## AIC = 56.4056, BIC = 69.7171 
    ## 
    ## Standard errors: MLE
    ## ------------------------------------------------------------------
    ##                     exp(Est.)     2.5%    97.5%    t val.        p
    ## ----------------- ----------- -------- -------- --------- --------
    ## (Intercept)            2.6436   0.8184   8.5395    1.6250   0.1057
    ## Gest_age               0.9407   0.9101   0.9723   -3.6234   0.0004
    ## FormulaTRUE            1.5304   1.0639   2.2014    2.2940   0.0228
    ## ------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.7258

``` r
summary(fit_meta2)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ Gest_age + Formula, family = Gamma(link = "log"), 
    ##     data = df_meta)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4290  -1.1854  -0.5520   0.1557   3.3715  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.97214    0.59826   1.625 0.105724    
    ## Gest_age    -0.06113    0.01687  -3.623 0.000367 ***
    ## FormulaTRUE  0.42553    0.18550   2.294 0.022813 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.725756)
    ## 
    ##     Null deviance: 369.30  on 205  degrees of freedom
    ## Residual deviance: 334.81  on 203  degrees of freedom
    ## AIC: 56.406
    ## 
    ## Number of Fisher Scoring iterations: 7

``` r
fit_meta3<-glm(SUM~SSU_counts+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta3, exp=TRUE, digits=4)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(3) = 41.2040, p = 0.0001
    ## Pseudo-R² (Cragg-Uhler) = 0.4544
    ## Pseudo-R² (McFadden) = 0.4109
    ## AIC = 53.2708, BIC = 69.9102 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------------
    ##                     exp(Est.)     2.5%     97.5%    t val.        p
    ## ----------------- ----------- -------- --------- --------- --------
    ## (Intercept)            3.2791   0.9505   11.3127    1.8796   0.0616
    ## SSU_counts             1.0000   1.0000    1.0000   -3.0192   0.0029
    ## Gest_age               0.9414   0.9091    0.9750   -3.3809   0.0009
    ## FormulaTRUE            1.6738   1.1388    2.4601    2.6215   0.0094
    ## -------------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.9173

``` r
summary(fit_meta3)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ SSU_counts + Gest_age + Formula, family = Gamma(link = "log"), 
    ##     data = df_meta)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4041  -1.1491  -0.5052   0.1877   3.5494  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.188e+00  6.318e-01   1.880 0.061602 .  
    ## SSU_counts  -8.269e-06  2.739e-06  -3.019 0.002861 ** 
    ## Gest_age    -6.035e-02  1.785e-02  -3.381 0.000867 ***
    ## FormulaTRUE  5.151e-01  1.965e-01   2.622 0.009420 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.917321)
    ## 
    ##     Null deviance: 369.30  on 205  degrees of freedom
    ## Residual deviance: 328.09  on 202  degrees of freedom
    ## AIC: 53.271
    ## 
    ## Number of Fisher Scoring iterations: 10

``` r
fit_meta4<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta4, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(4) = 61.92, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = 0.67
    ## Pseudo-R² (McFadden) = 0.63
    ## AIC = 38.86, BIC = 58.82 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------
    ##                     exp(Est.)   2.5%    97.5%   t val.      p
    ## ----------------- ----------- ------ -------- -------- ------
    ## (Intercept)             47.43   8.32   270.28     4.35   0.00
    ## SSU_counts               1.00   1.00     1.00    -4.10   0.00
    ## Age                      0.93   0.90     0.97    -3.99   0.00
    ## Gest_age                 0.89   0.85     0.93    -5.41   0.00
    ## FormulaTRUE              1.66   1.18     2.33     2.90   0.00
    ## -------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.51

``` r
summary(fit_meta4)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ SSU_counts + Age + Gest_age + Formula, family = Gamma(link = "log"), 
    ##     data = df_meta)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4640  -1.1051  -0.4923   0.2087   2.8770  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.859e+00  8.879e-01   4.346 2.20e-05 ***
    ## SSU_counts  -9.991e-06  2.434e-06  -4.104 5.90e-05 ***
    ## Age         -6.952e-02  1.741e-02  -3.993 9.15e-05 ***
    ## Gest_age    -1.204e-01  2.226e-02  -5.410 1.78e-07 ***
    ## FormulaTRUE  5.052e-01  1.743e-01   2.898  0.00417 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.507904)
    ## 
    ##     Null deviance: 369.30  on 205  degrees of freedom
    ## Residual deviance: 307.38  on 201  degrees of freedom
    ## AIC: 38.856
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
fit_meta4_AB<-glm(SUM~SSU_counts+Age+Gest_age+Formula, data=df_meta[df_meta$Antibiotics==TRUE,], family="Gamma"(link="log"))
summ(fit_meta4_AB, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 104
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(4) = 32.63, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = 0.40
    ## Pseudo-R² (McFadden) = 0.30
    ## AIC = 80.46, BIC = 96.32 
    ## 
    ## Standard errors: MLE
    ## -------------------------------------------------------------
    ##                     exp(Est.)   2.5%    97.5%   t val.      p
    ## ----------------- ----------- ------ -------- -------- ------
    ## (Intercept)             34.45   2.79   425.54     2.76   0.01
    ## SSU_counts               1.00   1.00     1.00    -2.86   0.01
    ## Age                      0.93   0.90     0.97    -3.43   0.00
    ## Gest_age                 0.89   0.82     0.97    -2.64   0.01
    ## FormulaTRUE              2.66   1.52     4.63     3.45   0.00
    ## -------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.43

``` r
summary(fit_meta4_AB)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ SSU_counts + Age + Gest_age + Formula, family = Gamma(link = "log"), 
    ##     data = df_meta[df_meta$Antibiotics == TRUE, ])
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.0118  -0.9711  -0.5388   0.2607   2.6646  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.539e+00  1.283e+00   2.759 0.006899 ** 
    ## SSU_counts  -1.436e-05  5.012e-06  -2.864 0.005102 ** 
    ## Age         -6.814e-02  1.985e-02  -3.433 0.000873 ***
    ## Gest_age    -1.154e-01  4.375e-02  -2.639 0.009670 ** 
    ## FormulaTRUE  9.774e-01  2.837e-01   3.445 0.000839 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.433909)
    ## 
    ##     Null deviance: 154.81  on 103  degrees of freedom
    ## Residual deviance: 122.17  on  99  degrees of freedom
    ## AIC: 80.458
    ## 
    ## Number of Fisher Scoring iterations: 11

``` r
fit_meta5<-glm(SUM~SSU_counts+Age+Gest_age+Antibiotics+Formula, data=df_meta, family="Gamma"(link="log"))
summ(fit_meta5, exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: SUM
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(5) = 61.92, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = 0.67
    ## Pseudo-R² (McFadden) = 0.63
    ## AIC = 40.85, BIC = 64.15 
    ## 
    ## Standard errors: MLE
    ## -----------------------------------------------------------------
    ##                         exp(Est.)   2.5%    97.5%   t val.      p
    ## --------------------- ----------- ------ -------- -------- ------
    ## (Intercept)                 49.54   4.02   609.96     3.05   0.00
    ## SSU_counts                   1.00   1.00     1.00    -4.01   0.00
    ## Age                          0.93   0.90     0.97    -3.91   0.00
    ## Gest_age                     0.89   0.83     0.94    -3.72   0.00
    ## AntibioticsTRUE              0.98   0.50     1.94    -0.05   0.96
    ## FormulaTRUE                  1.66   1.18     2.35     2.89   0.00
    ## -----------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 1.52

``` r
summary(fit_meta5)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ SSU_counts + Age + Gest_age + Antibiotics + 
    ##     Formula, family = Gamma(link = "log"), data = df_meta)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.4636  -1.1037  -0.4926   0.2052   2.8850  
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      3.903e+00  1.281e+00   3.047 0.002624 ** 
    ## SSU_counts      -9.974e-06  2.489e-06  -4.007 8.67e-05 ***
    ## Age             -6.943e-02  1.777e-02  -3.908 0.000128 ***
    ## Gest_age        -1.216e-01  3.271e-02  -3.717 0.000262 ***
    ## AntibioticsTRUE -1.656e-02  3.455e-01  -0.048 0.961833    
    ## FormulaTRUE      5.082e-01  1.761e-01   2.886 0.004336 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 1.516428)
    ## 
    ##     Null deviance: 369.30  on 205  degrees of freedom
    ## Residual deviance: 307.37  on 200  degrees of freedom
    ## AIC: 40.853
    ## 
    ## Number of Fisher Scoring iterations: 9

``` r
# Example chisquared test
anova(fit_meta4, fit_meta5, test="Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: SUM ~ SSU_counts + Age + Gest_age + Formula
    ## Model 2: SUM ~ SSU_counts + Age + Gest_age + Antibiotics + Formula
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1       201     307.38                     
    ## 2       200     307.37  1 0.004006    0.959

# Plot GLM

``` r
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
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
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

gest.age.plot_ff
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

# Meta-analysis Metaphlan2

``` r
otu_table_metaphlan<-read.table("merged_abundance_table_species_meta.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
tax_table_metaphlan<-read.table("tax_table_meta.txt", sep="\t", header=FALSE, fill = TRUE, row.names = 1, na.strings = "")
PHY_meta_metaphlan<-phyloseq(otu_table((otu_table_metaphlan), taxa_are_rows=TRUE), sample_data=sample_data(merged_phyloseq), tax_table(as.matrix(tax_table_metaphlan)))

colnames(tax_table(PHY_meta_metaphlan)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

toptaxa_gen<-find.top.taxa(PHY_meta_metaphlan,"Genus")
toptaxa_sp<-find.top.taxa(PHY_meta_metaphlan,"Species")
toptaxa_fam<-find.top.taxa(PHY_meta_metaphlan,"Family")


PHY_meta_metaphlan<-subset_taxa(PHY_meta_metaphlan, Domain=="k__Bacteria")

# Turn into relative data
PHY_meta_metaphlan<-transform_sample_counts(PHY_meta_metaphlan, function(x) x/sum(x))
sample_data(PHY_meta_metaphlan)$Top_gen<-toptaxa_gen$Genus
sample_data(PHY_meta_metaphlan)$Top_sp<-toptaxa_sp$Species
sample_data(PHY_meta_metaphlan)$Top_fam<-toptaxa_fam$Family


PHY_meta_metaphlan<-subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))


sample_data(PHY_meta_metaphlan)$Top_gen[!sample_data(PHY_meta_metaphlan)$Top_gen%in%c("g__Escherichia", "g__Klebsiella", "g__Bifidobacterium", "g__Bacteroides", "g__Staphylococcus", "g__Enterococcus", "g__Enterobacter", "g__Veillonella", "g__Haemophilus", "g__Clostridium")]<-"Other"


# Ordinate
temp<-sqrt(otu_table(PHY_meta_metaphlan))

merged_phyloseq_metaphlan_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaphlan)), tax_table(tax_table((PHY_meta_metaphlan))))

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
```

    ## Coordinate system already present. Adding new coordinate system, which will replace the existing one.

``` r
metaphlan.meta.ord.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# Meta-analysis Metaxa2

``` r
# Make Phyloseq object and modify tax table
otu_table_metaxa<-read.table("metaxa_meta-analysis_level6_otutab.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
tax_table_metaxa<-read.table("metaxa_meta-analysis_level6_tax.txt", sep="\t", header=FALSE, fill = TRUE, row.names = 1, na.strings = "")
tax_table_metaxa[is.na(tax_table_metaxa)]<-"Unclassified"

tax_table_metaxa$V8<-paste(tax_table_metaxa$V3, tax_table_metaxa$V4, tax_table_metaxa$V5,  tax_table_metaxa$V6, sep=";")
tax_table_metaxa$V9<-paste(tax_table_metaxa$V3, tax_table_metaxa$V4, tax_table_metaxa$V5,  tax_table_metaxa$V6,  tax_table_metaxa$V7, sep=";")



PHY_meta_metaxa<-phyloseq(otu_table((otu_table_metaxa), taxa_are_rows=TRUE), sample_data=sample_data(merged_phyloseq), tax_table(as.matrix(tax_table_metaxa)))

# Change column names to somethig meaningful
colnames(tax_table(PHY_meta_metaxa)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Fam_path", "Gen_path")

PHY_meta_metaxa<-subset_taxa(PHY_meta_metaxa, Domain=="Bacteria")
# Turn into relative data
PHY_meta_metaxa_rel<-transform_sample_counts(PHY_meta_metaxa, function(x) x/sum(x))

# Add dominant taxa
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

#Modify sample data to include dominant taxa
sample_data(PHY_meta_metaxa_rel)$Top_fam_meta<-toptaxa$Family
sample_data(PHY_meta_metaxa_rel)$Top_gen_meta<-paste(toptaxa_gen$Genus, "-" ,toptaxa$Family)
sample_data(PHY_meta_metaxa_rel)$Top_gen_meta[!sample_data(PHY_meta_metaxa_rel)$Top_gen_meta%in%c("Unclassified - Enterobacteriaceae", "Bacteroides - Bacteroidaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Bifidobacterium - Bifidobacteriaceae", "Enterococcus - Enterococcacea", "Clostridium - Clostridiaceae",  "Streptococcus - Streptococcaceae", "Veillonella - Veillonellaceae")]<-"Other"
sample_data(PHY_meta_metaxa_rel)$Top_fam_meta[!sample_data(PHY_meta_metaxa_rel)$Top_fam_meta%in%c("Enterobacteriaceae", "Bacteroidaceae", "Staphylococcaceae", "Bifidobacteriaceae", "Enterococcaceae", "Clostridiaceae")]<-"Other"
```

# Ordination

``` r
# Ordinate
# Metaxa
# Squareroot transformation

temp<-sqrt(otu_table(PHY_meta_metaxa_rel))
merged_phyloseq_metaxa_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaxa_rel)), tax_table(tax_table((PHY_meta_metaxa_rel))))

# Ordinate

PHY_SP_meta_ord<-ordinate(merged_phyloseq_metaxa_sqrt, method="PCoA", distance="horn")

# Plot

p<-plot_ordination(merged_phyloseq_metaxa_sqrt, PHY_SP_meta_ord, color="Top_gen_meta")
p$layers <- p$layers[-1]

metaxa.meta.ord.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", guide=guide_legend(nrow=3),breaks=c("Bacteroides - Bacteroidaceae", "Bifidobacterium - Bifidobacteriaceae", "Clostridium - Clostridiaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Streptococcus - Streptococcaceae", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae", "Other"), labels=c("Bacteroides", "Bifidobacterium", "Clostridium", "Escherichia-Shigella", "Staphylococcus", "Streptococcus", "Uncl. Enterobacteriaceae", "Veillonella", "Other")) + geom_point(size=4)  +theme_classic() + labs(title="Microbial community") + theme(text = element_text(size=16), legend.position = "bottom") +
  stat_ellipse(level = 0.90, linetype=2, size=1.5,show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
   ylim(-0.75, 0.75) +
  xlim(-0.75, 0.75)


# ARGs
# Squareroot transformation
temp<-sqrt(otu_table(merged_phyloseq))
merged_phyloseq_sqrt<-phyloseq(otu_table(temp), sample_data(sample_data(PHY_meta_metaxa_rel)), tax_table(tax_table((merged_phyloseq))))
# Ordinate
PHY_ARG_meta_ord<-ordinate(merged_phyloseq_sqrt, method="PCoA", distance="horn")

# Plot
p<-plot_ordination(merged_phyloseq_sqrt, PHY_ARG_meta_ord, color="Top_gen_meta")
p$layers <- p$layers[-1]

ARG.meta.ord.plot<-p+scale_color_brewer(palette="BrBG", "Dominant genus", breaks=c("Bacteroides - Bacteroidaceae", "Bifidobacterium - Bifidobacteriaceae", "Clostridium - Clostridiaceae", "Escherichia-Shigella - Enterobacteriaceae", "Staphylococcus - Staphylococcaceae", "Streptococcus - Streptococcaceae", "Unclassified - Enterobacteriaceae", "Veillonella - Veillonellaceae", "Other")) + geom_point(size=4)  +theme_classic() + labs(title="ARGs") + theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "pt")) +
  stat_ellipse(level = 0.90, size=1.5, linetype=2, show.legend = FALSE) +
  coord_equal() +
  geom_point(shape=21, size=4, color="black") +
   ylim(-0.75, 0.75) +
  xlim(-0.75, 0.75) #+
  coord_equal() 
```

    ## <ggproto object: Class CoordFixed, CoordCartesian, Coord, gg>
    ##     aspect: function
    ##     backtransform_range: function
    ##     clip: on
    ##     default: FALSE
    ##     distance: function
    ##     expand: TRUE
    ##     is_free: function
    ##     is_linear: function
    ##     labels: function
    ##     limits: list
    ##     modify_scales: function
    ##     range: function
    ##     ratio: 1
    ##     render_axis_h: function
    ##     render_axis_v: function
    ##     render_bg: function
    ##     render_fg: function
    ##     setup_data: function
    ##     setup_layout: function
    ##     setup_panel_params: function
    ##     setup_params: function
    ##     transform: function
    ##     super:  <ggproto object: Class CoordFixed, CoordCartesian, Coord, gg>

``` r
# Test differences between infants dominated by different genera
  
df_temp<-data.frame(SUM=sample_sums(merged_phyloseq), Top_genus=sample_data(merged_phyloseq_sqrt)$Top_gen_meta, Formula=sample_data(merged_phyloseq_sqrt)$Formula2)

fit_meta_genera<-glm(SUM~Top_genus, data=df_temp, family="Gamma"(link="log"))

# Post hoc test for the model

glht.mod <- glht(fit_meta_genera, mcp(Top_genus = "Tukey"))
summary(glht(glht.mod))
```

    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored
    
    ## Warning in chkdots(...): Argument(s) 'complete' passed to '...' are ignored

    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Linear Hypotheses:
    ##                                                                                       Estimate
    ## Bifidobacterium - Bifidobacteriaceae - Bacteroides - Bacteroidaceae == 0              -0.54796
    ## Clostridium - Clostridiaceae - Bacteroides - Bacteroidaceae == 0                      -0.96819
    ## Escherichia-Shigella - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0          0.76813
    ## Other - Bacteroides - Bacteroidaceae == 0                                             -0.20258
    ## Staphylococcus - Staphylococcaceae - Bacteroides - Bacteroidaceae == 0                 1.20732
    ## Streptococcus - Streptococcaceae - Bacteroides - Bacteroidaceae == 0                   0.11017
    ## Unclassified - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0                 -0.17905
    ## Veillonella - Veillonellaceae - Bacteroides - Bacteroidaceae == 0                     -0.86043
    ## Clostridium - Clostridiaceae - Bifidobacterium - Bifidobacteriaceae == 0              -0.42023
    ## Escherichia-Shigella - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0  1.31609
    ## Other - Bifidobacterium - Bifidobacteriaceae == 0                                      0.34538
    ## Staphylococcus - Staphylococcaceae - Bifidobacterium - Bifidobacteriaceae == 0         1.75527
    ## Streptococcus - Streptococcaceae - Bifidobacterium - Bifidobacteriaceae == 0           0.65813
    ## Unclassified - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0          0.36890
    ## Veillonella - Veillonellaceae - Bifidobacterium - Bifidobacteriaceae == 0             -0.31247
    ## Escherichia-Shigella - Enterobacteriaceae - Clostridium - Clostridiaceae == 0          1.73632
    ## Other - Clostridium - Clostridiaceae == 0                                              0.76560
    ## Staphylococcus - Staphylococcaceae - Clostridium - Clostridiaceae == 0                 2.17550
    ## Streptococcus - Streptococcaceae - Clostridium - Clostridiaceae == 0                   1.07835
    ## Unclassified - Enterobacteriaceae - Clostridium - Clostridiaceae == 0                  0.78913
    ## Veillonella - Veillonellaceae - Clostridium - Clostridiaceae == 0                      0.10776
    ## Other - Escherichia-Shigella - Enterobacteriaceae == 0                                -0.97071
    ## Staphylococcus - Staphylococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0    0.43919
    ## Streptococcus - Streptococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0     -0.65796
    ## Unclassified - Enterobacteriaceae - Escherichia-Shigella - Enterobacteriaceae == 0    -0.94718
    ## Veillonella - Veillonellaceae - Escherichia-Shigella - Enterobacteriaceae == 0        -1.62856
    ## Staphylococcus - Staphylococcaceae - Other == 0                                        1.40990
    ## Streptococcus - Streptococcaceae - Other == 0                                          0.31275
    ## Unclassified - Enterobacteriaceae - Other == 0                                         0.02353
    ## Veillonella - Veillonellaceae - Other == 0                                            -0.65785
    ## Streptococcus - Streptococcaceae - Staphylococcus - Staphylococcaceae == 0            -1.09715
    ## Unclassified - Enterobacteriaceae - Staphylococcus - Staphylococcaceae == 0           -1.38637
    ## Veillonella - Veillonellaceae - Staphylococcus - Staphylococcaceae == 0               -2.06775
    ## Unclassified - Enterobacteriaceae - Streptococcus - Streptococcaceae == 0             -0.28922
    ## Veillonella - Veillonellaceae - Streptococcus - Streptococcaceae == 0                 -0.97060
    ## Veillonella - Veillonellaceae - Unclassified - Enterobacteriaceae == 0                -0.68138
    ##                                                                                       Std. Error
    ## Bifidobacterium - Bifidobacteriaceae - Bacteroides - Bacteroidaceae == 0                 0.37903
    ## Clostridium - Clostridiaceae - Bacteroides - Bacteroidaceae == 0                         0.40894
    ## Escherichia-Shigella - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0            0.29722
    ## Other - Bacteroides - Bacteroidaceae == 0                                                0.25325
    ## Staphylococcus - Staphylococcaceae - Bacteroides - Bacteroidaceae == 0                   0.30524
    ## Streptococcus - Streptococcaceae - Bacteroides - Bacteroidaceae == 0                     0.45195
    ## Unclassified - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0                    0.22460
    ## Veillonella - Veillonellaceae - Bacteroides - Bacteroidaceae == 0                        0.52049
    ## Clostridium - Clostridiaceae - Bifidobacterium - Bifidobacteriaceae == 0                 0.48543
    ## Escherichia-Shigella - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0    0.39592
    ## Other - Bifidobacterium - Bifidobacteriaceae == 0                                        0.36407
    ## Staphylococcus - Staphylococcaceae - Bifidobacterium - Bifidobacteriaceae == 0           0.40197
    ## Streptococcus - Streptococcaceae - Bifidobacterium - Bifidobacteriaceae == 0             0.52218
    ## Unclassified - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0            0.34475
    ## Veillonella - Veillonellaceae - Bifidobacterium - Bifidobacteriaceae == 0                0.58251
    ## Escherichia-Shigella - Enterobacteriaceae - Clostridium - Clostridiaceae == 0            0.42464
    ## Other - Clostridium - Clostridiaceae == 0                                                0.39511
    ## Staphylococcus - Staphylococcaceae - Clostridium - Clostridiaceae == 0                   0.43029
    ## Streptococcus - Streptococcaceae - Clostridium - Clostridiaceae == 0                     0.54427
    ## Unclassified - Enterobacteriaceae - Clostridium - Clostridiaceae == 0                    0.37739
    ## Veillonella - Veillonellaceae - Clostridium - Clostridiaceae == 0                        0.60240
    ## Other - Escherichia-Shigella - Enterobacteriaceae == 0                                   0.27789
    ## Staphylococcus - Staphylococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0      0.32597
    ## Streptococcus - Streptococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0        0.46620
    ## Unclassified - Enterobacteriaceae - Escherichia-Shigella - Enterobacteriaceae == 0       0.25205
    ## Veillonella - Veillonellaceae - Escherichia-Shigella - Enterobacteriaceae == 0           0.53291
    ## Staphylococcus - Staphylococcaceae - Other == 0                                          0.28645
    ## Streptococcus - Streptococcaceae - Other == 0                                            0.43948
    ## Unclassified - Enterobacteriaceae - Other == 0                                           0.19831
    ## Veillonella - Veillonellaceae - Other == 0                                               0.50970
    ## Streptococcus - Streptococcaceae - Staphylococcus - Staphylococcaceae == 0               0.47135
    ## Unclassified - Enterobacteriaceae - Staphylococcus - Staphylococcaceae == 0              0.26146
    ## Veillonella - Veillonellaceae - Staphylococcus - Staphylococcaceae == 0                  0.53743
    ## Unclassified - Enterobacteriaceae - Streptococcus - Streptococcaceae == 0                0.42361
    ## Veillonella - Veillonellaceae - Streptococcus - Streptococcaceae == 0                    0.63239
    ## Veillonella - Veillonellaceae - Unclassified - Enterobacteriaceae == 0                   0.49609
    ##                                                                                       z value
    ## Bifidobacterium - Bifidobacteriaceae - Bacteroides - Bacteroidaceae == 0               -1.446
    ## Clostridium - Clostridiaceae - Bacteroides - Bacteroidaceae == 0                       -2.368
    ## Escherichia-Shigella - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0           2.584
    ## Other - Bacteroides - Bacteroidaceae == 0                                              -0.800
    ## Staphylococcus - Staphylococcaceae - Bacteroides - Bacteroidaceae == 0                  3.955
    ## Streptococcus - Streptococcaceae - Bacteroides - Bacteroidaceae == 0                    0.244
    ## Unclassified - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0                  -0.797
    ## Veillonella - Veillonellaceae - Bacteroides - Bacteroidaceae == 0                      -1.653
    ## Clostridium - Clostridiaceae - Bifidobacterium - Bifidobacteriaceae == 0               -0.866
    ## Escherichia-Shigella - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0   3.324
    ## Other - Bifidobacterium - Bifidobacteriaceae == 0                                       0.949
    ## Staphylococcus - Staphylococcaceae - Bifidobacterium - Bifidobacteriaceae == 0          4.367
    ## Streptococcus - Streptococcaceae - Bifidobacterium - Bifidobacteriaceae == 0            1.260
    ## Unclassified - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0           1.070
    ## Veillonella - Veillonellaceae - Bifidobacterium - Bifidobacteriaceae == 0              -0.536
    ## Escherichia-Shigella - Enterobacteriaceae - Clostridium - Clostridiaceae == 0           4.089
    ## Other - Clostridium - Clostridiaceae == 0                                               1.938
    ## Staphylococcus - Staphylococcaceae - Clostridium - Clostridiaceae == 0                  5.056
    ## Streptococcus - Streptococcaceae - Clostridium - Clostridiaceae == 0                    1.981
    ## Unclassified - Enterobacteriaceae - Clostridium - Clostridiaceae == 0                   2.091
    ## Veillonella - Veillonellaceae - Clostridium - Clostridiaceae == 0                       0.179
    ## Other - Escherichia-Shigella - Enterobacteriaceae == 0                                 -3.493
    ## Staphylococcus - Staphylococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0     1.347
    ## Streptococcus - Streptococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0      -1.411
    ## Unclassified - Enterobacteriaceae - Escherichia-Shigella - Enterobacteriaceae == 0     -3.758
    ## Veillonella - Veillonellaceae - Escherichia-Shigella - Enterobacteriaceae == 0         -3.056
    ## Staphylococcus - Staphylococcaceae - Other == 0                                         4.922
    ## Streptococcus - Streptococcaceae - Other == 0                                           0.712
    ## Unclassified - Enterobacteriaceae - Other == 0                                          0.119
    ## Veillonella - Veillonellaceae - Other == 0                                             -1.291
    ## Streptococcus - Streptococcaceae - Staphylococcus - Staphylococcaceae == 0             -2.328
    ## Unclassified - Enterobacteriaceae - Staphylococcus - Staphylococcaceae == 0            -5.302
    ## Veillonella - Veillonellaceae - Staphylococcus - Staphylococcaceae == 0                -3.847
    ## Unclassified - Enterobacteriaceae - Streptococcus - Streptococcaceae == 0              -0.683
    ## Veillonella - Veillonellaceae - Streptococcus - Streptococcaceae == 0                  -1.535
    ## Veillonella - Veillonellaceae - Unclassified - Enterobacteriaceae == 0                 -1.374
    ##                                                                                       Pr(>|z|)
    ## Bifidobacterium - Bifidobacteriaceae - Bacteroides - Bacteroidaceae == 0               0.86271
    ## Clostridium - Clostridiaceae - Bacteroides - Bacteroidaceae == 0                       0.27298
    ## Escherichia-Shigella - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0          0.17173
    ## Other - Bacteroides - Bacteroidaceae == 0                                              0.99619
    ## Staphylococcus - Staphylococcaceae - Bacteroides - Bacteroidaceae == 0                 0.00211
    ## Streptococcus - Streptococcaceae - Bacteroides - Bacteroidaceae == 0                   1.00000
    ## Unclassified - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0                  0.99628
    ## Veillonella - Veillonellaceae - Bacteroides - Bacteroidaceae == 0                      0.74896
    ## Clostridium - Clostridiaceae - Bifidobacterium - Bifidobacteriaceae == 0               0.99345
    ## Escherichia-Shigella - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0  0.02144
    ## Other - Bifidobacterium - Bifidobacteriaceae == 0                                      0.98796
    ## Staphylococcus - Staphylococcaceae - Bifidobacterium - Bifidobacteriaceae == 0         < 0.001
    ## Streptococcus - Streptococcaceae - Bifidobacterium - Bifidobacteriaceae == 0           0.93281
    ## Unclassified - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0          0.97418
    ## Veillonella - Veillonellaceae - Bifidobacterium - Bifidobacteriaceae == 0              0.99979
    ## Escherichia-Shigella - Enterobacteriaceae - Clostridium - Clostridiaceae == 0          0.00116
    ## Other - Clostridium - Clostridiaceae == 0                                              0.55283
    ## Staphylococcus - Staphylococcaceae - Clostridium - Clostridiaceae == 0                 < 0.001
    ## Streptococcus - Streptococcaceae - Clostridium - Clostridiaceae == 0                   0.52198
    ## Unclassified - Enterobacteriaceae - Clostridium - Clostridiaceae == 0                  0.44531
    ## Veillonella - Veillonellaceae - Clostridium - Clostridiaceae == 0                      1.00000
    ## Other - Escherichia-Shigella - Enterobacteriaceae == 0                                 0.01213
    ## Staphylococcus - Staphylococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0    0.90380
    ## Streptococcus - Streptococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0      0.87791
    ## Unclassified - Enterobacteriaceae - Escherichia-Shigella - Enterobacteriaceae == 0     0.00438
    ## Veillonella - Veillonellaceae - Escherichia-Shigella - Enterobacteriaceae == 0         0.04886
    ## Staphylococcus - Staphylococcaceae - Other == 0                                        < 0.001
    ## Streptococcus - Streptococcaceae - Other == 0                                          0.99833
    ## Unclassified - Enterobacteriaceae - Other == 0                                         1.00000
    ## Veillonella - Veillonellaceae - Other == 0                                             0.92336
    ## Streptococcus - Streptococcaceae - Staphylococcus - Staphylococcaceae == 0             0.29588
    ## Unclassified - Enterobacteriaceae - Staphylococcus - Staphylococcaceae == 0            < 0.001
    ## Veillonella - Veillonellaceae - Staphylococcus - Staphylococcaceae == 0                0.00330
    ## Unclassified - Enterobacteriaceae - Streptococcus - Streptococcaceae == 0              0.99876
    ## Veillonella - Veillonellaceae - Streptococcus - Streptococcaceae == 0                  0.81789
    ## Veillonella - Veillonellaceae - Unclassified - Enterobacteriaceae == 0                 0.89382
    ##                                                                                          
    ## Bifidobacterium - Bifidobacteriaceae - Bacteroides - Bacteroidaceae == 0                 
    ## Clostridium - Clostridiaceae - Bacteroides - Bacteroidaceae == 0                         
    ## Escherichia-Shigella - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0            
    ## Other - Bacteroides - Bacteroidaceae == 0                                                
    ## Staphylococcus - Staphylococcaceae - Bacteroides - Bacteroidaceae == 0                ** 
    ## Streptococcus - Streptococcaceae - Bacteroides - Bacteroidaceae == 0                     
    ## Unclassified - Enterobacteriaceae - Bacteroides - Bacteroidaceae == 0                    
    ## Veillonella - Veillonellaceae - Bacteroides - Bacteroidaceae == 0                        
    ## Clostridium - Clostridiaceae - Bifidobacterium - Bifidobacteriaceae == 0                 
    ## Escherichia-Shigella - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0 *  
    ## Other - Bifidobacterium - Bifidobacteriaceae == 0                                        
    ## Staphylococcus - Staphylococcaceae - Bifidobacterium - Bifidobacteriaceae == 0        ***
    ## Streptococcus - Streptococcaceae - Bifidobacterium - Bifidobacteriaceae == 0             
    ## Unclassified - Enterobacteriaceae - Bifidobacterium - Bifidobacteriaceae == 0            
    ## Veillonella - Veillonellaceae - Bifidobacterium - Bifidobacteriaceae == 0                
    ## Escherichia-Shigella - Enterobacteriaceae - Clostridium - Clostridiaceae == 0         ** 
    ## Other - Clostridium - Clostridiaceae == 0                                                
    ## Staphylococcus - Staphylococcaceae - Clostridium - Clostridiaceae == 0                ***
    ## Streptococcus - Streptococcaceae - Clostridium - Clostridiaceae == 0                     
    ## Unclassified - Enterobacteriaceae - Clostridium - Clostridiaceae == 0                    
    ## Veillonella - Veillonellaceae - Clostridium - Clostridiaceae == 0                        
    ## Other - Escherichia-Shigella - Enterobacteriaceae == 0                                *  
    ## Staphylococcus - Staphylococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0      
    ## Streptococcus - Streptococcaceae - Escherichia-Shigella - Enterobacteriaceae == 0        
    ## Unclassified - Enterobacteriaceae - Escherichia-Shigella - Enterobacteriaceae == 0    ** 
    ## Veillonella - Veillonellaceae - Escherichia-Shigella - Enterobacteriaceae == 0        *  
    ## Staphylococcus - Staphylococcaceae - Other == 0                                       ***
    ## Streptococcus - Streptococcaceae - Other == 0                                            
    ## Unclassified - Enterobacteriaceae - Other == 0                                           
    ## Veillonella - Veillonellaceae - Other == 0                                               
    ## Streptococcus - Streptococcaceae - Staphylococcus - Staphylococcaceae == 0               
    ## Unclassified - Enterobacteriaceae - Staphylococcus - Staphylococcaceae == 0           ***
    ## Veillonella - Veillonellaceae - Staphylococcus - Staphylococcaceae == 0               ** 
    ## Unclassified - Enterobacteriaceae - Streptococcus - Streptococcaceae == 0                
    ## Veillonella - Veillonellaceae - Streptococcus - Streptococcaceae == 0                    
    ## Veillonella - Veillonellaceae - Unclassified - Enterobacteriaceae == 0                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
# Adonis
# ARGS
arg_adonis_meta<-pairwise.adonis(x=(t(otu_table(merged_phyloseq_sqrt))), factors = sample_data(merged_phyloseq_sqrt)$Formula2,  sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

``` r
arg_adonis_meta
```

    ##           pairs F.Model         R2 p.value p.adjusted sig
    ## 1 TRUE vs FALSE 2.68494 0.01106348   0.002      0.002   *

``` r
# Adonis
# Species
species_adonis_meta<-pairwise.adonis(x=(t(otu_table(merged_phyloseq_metaphlan_sqrt))), factors = sample_data(merged_phyloseq_metaphlan_sqrt)$Formula2,  sim.function = "vegdist", sim.method = "horn", p.adjust.m = "fdr")
```

    ## [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

``` r
species_adonis_meta
```

    ##           pairs  F.Model         R2 p.value p.adjusted sig
    ## 1 TRUE vs FALSE 3.675967 0.01770049   0.002      0.002   *

# DESeq meta-analysis ARGs

``` r
# ARGs and breastfeeding exclusively vs formula
temp_arg<-otu_table(subset_samples(merged_phyloseq, !(Formula2==FALSE&Fortifier==TRUE)))*5*10^4+1
ARG_tax2<-ARG_tax
# Merge TEM genes
ARG_tax2$V3<-gsub('TEM-[[:digit:]]+', 'TEM', ARG_tax2$V3)

# Modify gene names and taxonomy table
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

# Create phyloseq object for DESeq
ARG_DSQ<-phyloseq(otu_table(temp_arg, taxa_are_rows=T), sample_data(sample_data(merged_phyloseq)), tax_table(as.matrix(ARG_tax2)))
ARG_DSQ_glom<-tax_glom(ARG_DSQ, taxrank="V3")
dds_arg = phyloseq_to_deseq2(ARG_DSQ_glom, ~ Formula2)
```

    ## converting counts to integer mode

``` r
dds_arg = DESeq(dds_arg, fitType = "mean", test ="Wald")
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 489 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
resultsNames(dds_arg)
```

    ## [1] "Intercept"    "Formula2TRUE"

``` r
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


# Plot
Deseq_meta<-ggplot(sigtab_arg, aes(x=V3, y=log2FoldChange, color=V2, fill=V2, size=100*baseMean)) + geom_jitter(alpha=0.8, width = 0.1, shape=21, color="black")+theme_minimal()+ theme(axis.title.y = element_text(size=rel(1.3)), legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)), panel.grid.major = element_blank(), legend.position="bottom", axis.text.x=element_blank())+ geom_text_repel(aes(label=V3), size=4, color= "black", fontface = "italic") +  scale_fill_brewer(palette="BrBG", "ARG class") + scale_size(range=c(3,20),breaks=c(3, 5, 10, 15, 20)) + labs(color="") +guides(size=FALSE) +  ylim(-6, 6) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") + xlab("") +
guides(fill = guide_legend(override.aes = list(size = 4)),
            size = guide_legend(override.aes = list(linetype = 0))) +
  guides(size = guide_legend(override.aes = list(linetype = 0)))

Deseq_meta
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

# DESeq meta-analysis Metaxa

``` r
# Metaxa
# Taxa and breastfeeding exclusively vs formula
# Pseudocounts
temp_mtx<-otu_table(PHY_meta_metaxa)+1
mtx_tax2<-tax_table(PHY_meta_metaxa)

mtx_DSQ<-phyloseq(otu_table(temp_mtx, taxa_are_rows=T), sample_data(PHY_meta_metaxa), tax_table(as.matrix(mtx_tax2)))
mtx_DSQ_glom<-tax_glom(mtx_DSQ, taxrank="Genus")
dds_mtx = phyloseq_to_deseq2(mtx_DSQ_glom, ~ Formula2)
```

    ## converting counts to integer mode

``` r
dds_mtx = DESeq(dds_mtx, fitType = "mean", test ="Wald")
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 110 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
res_mtx = results(dds_mtx, cooksCutoff = FALSE)
alpha = 0.05
sigtab_mtx = res_mtx[which(res_mtx$padj < alpha), ]
sigtab_mtx = cbind(as(sigtab_mtx, "data.frame"), as(tax_table(mtx_DSQ)[rownames(sigtab_mtx), ], "matrix"))

sigtab_mtx<-sigtab_mtx[1:12]
sigtab_mtx<-sigtab_mtx[sigtab_mtx$baseMean>10,]
sigtab_mtx<-sigtab_mtx[!sigtab_mtx$Family%in%c("Unclassified", "Family"),]

# Plot

pal<-wes_palette("IsleofDogs1", 6, type="discrete")

# Plot
Deseq_meta_mtx<-ggplot(sigtab_mtx, aes(x=Family, y=log2FoldChange, color=Class, fill=Class, size=100*baseMean)) + geom_jitter(alpha=0.9, width = 0.2, shape=21, color="black")+theme_minimal()+ theme(axis.text.x = element_text(angle =45, hjust = 0.9, vjust = 0.97, size=rel(1.3), face="italic"), axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_text(size=rel(1.3)), legend.text = element_text(size=rel(1.2)), legend.title = element_text(size=rel(1.2)), legend.position="bottom")+ geom_text_repel(aes(label=Genus, fontface="italic"), size=4, color= "black") +  scale_fill_manual(values = pal) + scale_size(range=c(3,20),breaks=c(3, 5, 10, 15, 20)) + labs(color="") +guides(size=FALSE) +  ylim(-5, 5) + geom_hline(yintercept = 0, linetype=2, color="grey", alpha=0.5) +ylab("log2 fold change\nBreast milk     Formula") +
guides(fill = guide_legend(override.aes = list(size = 4)),
            size = guide_legend(override.aes = list(linetype = 0)))

Deseq_meta_mtx
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

# Differences in ARG sums related to top genus using GLMs

``` r
# Check differences in ARG sums related to top genus using GLMs
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

    ## Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position

``` r
ARG.gen.sum.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

# Diversity

``` r
# Calculate diversity and make dataframe
df <- data.frame(DIV=diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "shannon"), DIV2=diversity(t(otu_table(subset_samples(PHY_meta_metaphlan))), index = "simpson"), Formula=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Formula2, Gest_age=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Gest_age,
                 Age=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Age,
                  Antibiotics=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Antibiotics,
                  Delivery_mode=sample_data(subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE)))$Delivery_mode)


summ(glm(DIV2~Age+Gest_age+Antibiotics+Formula, data=df, family=Gamma(link="log")), exp=TRUE)
```

    ## MODEL INFO:
    ## Observations: 206
    ## Dependent Variable: DIV2
    ## Type: Generalized linear model
    ##   Family: Gamma 
    ##   Link function: log 
    ## 
    ## MODEL FIT:
    ## χ²(4) = 7.11, p = 0.00
    ## Pseudo-R² (Cragg-Uhler) = 0.15
    ## Pseudo-R² (McFadden) = 0.13
    ## AIC = 83.47, BIC = 103.43 
    ## 
    ## Standard errors: MLE
    ## ----------------------------------------------------------------
    ##                         exp(Est.)   2.5%   97.5%   t val.      p
    ## --------------------- ----------- ------ ------- -------- ------
    ## (Intercept)                  0.14   0.05    0.42    -3.48   0.00
    ## Age                          1.03   1.01    1.04     3.43   0.00
    ## Gest_age                     1.03   1.00    1.06     2.22   0.03
    ## AntibioticsTRUE              1.01   0.75    1.35     0.04   0.97
    ## FormulaTRUE                  0.78   0.67    0.91    -3.13   0.00
    ## ----------------------------------------------------------------
    ## 
    ## Estimated dispersion parameter = 0.3

``` r
df$Formula[df$Formula==TRUE]<-"Yes"
df$Formula[df$Formula==FALSE]<-"No"

# Plot 
metaphlan.shandiv.plot<-ggplot(df, aes(x=Formula, y=DIV, fill=Formula))  +geom_jitter(width=0.4, size=4, shape=21) + geom_boxplot(alpha=0.5) + scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes"))+ theme_classic()  +guides(alpha=FALSE, fill=FALSE, color=FALSE)+
 labs(y="Shannon diversity")  + ggtitle("Shannon diversity of species") 

metaphlan.simdiv.plot<-ggplot(df, aes(x=Formula, y=DIV2, fill=Formula))  +geom_jitter(width=0.4, size=4, shape=21) + geom_boxplot(alpha=0.5) + scale_fill_manual(values=c("darkgoldenrod3", "azure1"), labels=c("No", "Yes"))+ theme_classic()  +guides(alpha=FALSE, fill=FALSE, color=FALSE)+
 labs(y="Simpson diversity")  + ggtitle("Simpson diversity of species") 

metaphlan.simdiv.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
meta.div.plot<-cowplot::plot_grid(metaphlan.shandiv.plot, metaphlan.simdiv.plot, labels = c("c", "d"))

meta.div.plot
```

![](Supplementary_software_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
# Analysis of variance

a0<-aov(diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "shannon")~as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data=as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
TukeyHSD(a0)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2 == FALSE & Fortifier == TRUE))))), index = "shannon") ~ as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data = as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
    ## 
    ## $`as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2)`
    ##                  diff        lwr        upr     p adj
    ## TRUE-FALSE -0.2636259 -0.4160732 -0.1111785 0.0007841

``` r
a1<-aov(diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2==FALSE&Fortifier==TRUE))))), index = "simpson")~as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data=as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
TukeyHSD(a1)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = diversity(t(otu_table((subset_samples(PHY_meta_metaphlan, !(Formula2 == FALSE & Fortifier == TRUE))))), index = "simpson") ~ as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2), data = as.data.frame((otu_table(subset_samples(PHY_meta_metaphlan)))))
    ## 
    ## $`as.factor(sample_data((subset_samples(PHY_meta_metaphlan)))$Formula2)`
    ##                  diff        lwr         upr    p adj
    ## TRUE-FALSE -0.1136569 -0.1797068 -0.04760699 0.000831

# Effect of donor milk

``` r
# Make phyloseq object for analysis of the effect of donor using cohorts that include donor milk
donor_phyloseq<-subset_samples(merged_phyloseq, Study%in%c("NEC", "Gibson"))

sample_data(donor_phyloseq)$Donor_Milk<-sample_data(donor_phyloseq)$Donor_Milk>0

# Make dataframe
df_don<-data.frame(SUM=sample_sums(donor_phyloseq),
                   Formula=sample_data(donor_phyloseq)$Formula,
                   Gest_age=sample_data(donor_phyloseq)$Gest_age,
                   Fortifier=sample_data(donor_phyloseq)$Fortifier,
                   Age=sample_data(donor_phyloseq)$Age,
                   Antibiotics=sample_data(donor_phyloseq)$Antibiotics,
                   Donor=sample_data(donor_phyloseq)$Donor_Milk,
                   Study=sample_data(donor_phyloseq)$Study)
                   
df_don<-df_don[df_don$Formula=="FALSE",]

# Fit GLM
fitD<-glm(SUM~Gest_age+Antibiotics+Donor, data=df_don, family=Gamma(link="log"))

summary(fitD)
```

    ## 
    ## Call:
    ## glm(formula = SUM ~ Gest_age + Antibiotics + Donor, family = Gamma(link = "log"), 
    ##     data = df_don)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8679  -0.7518  -0.3046   0.2915   2.5529  
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      3.41265    1.34274   2.542 0.013509 *  
    ## Gest_age        -0.16295    0.04277  -3.810 0.000318 ***
    ## AntibioticsTRUE  0.09491    0.40620   0.234 0.816002    
    ## DonorTRUE       -0.03821    0.35259  -0.108 0.914058    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.9231876)
    ## 
    ##     Null deviance: 65.627  on 66  degrees of freedom
    ## Residual deviance: 49.253  on 63  degrees of freedom
    ## AIC: -27.08
    ## 
    ## Number of Fisher Scoring iterations: 7
