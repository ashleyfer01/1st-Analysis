#Analysis of the Differntial Gene expression between TB+AMs and TB+IMs 
library(tidyverse) 
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
raw<-read.csv("lung_rawdata.csv",header = T) 
genenames=raw[,1]
mymatrix_numeric<-raw[1:23361,2:13]#numeric matrix of all the gene counts 
rownames(mymatrix_numeric)<-c(genenames)#adding the list of gene names to my numeric matrix 
raw_DGEList<-DGEList(counts =mymatrix_numeric)
raw_DGEList
cpm_from_rawDGElist<-cpm(raw_DGEList)
colSums(cpm_from_rawDGElist)
#study design matrix 
sample_id<-c("AM-TB1","AM-TB2","AM-TB3","AM+TB1","AM+TB2","AM+TB3","IM-TB1","IM-TB2","IM-TB3","IM+TB1","IM+TB2","IM+TB3")
cell_type<-c("AM","AM","AM","AM","AM","AM","IM","IM","IM","IM","IM","IM")
group<-c("AM_TB_Neg","AM_TB_Neg","AM_TB_Neg","AM_TB_Pos","AM_TB_Pos","AM_TB_Pos","IM_TB_Neg","IM_TB_Neg","IM_TB_Neg","IM_TB_Pos","IM_TB_Pos","IM_TB_Pos")
study_design<-data.frame(sample_id,cell_type,group)
study_design
sample_names2<-c("AMsample1","AMsample2","AMsample3","AMsample4","AMsample5","AMsample6","IMsample1","IMsample2","IMsample3","IMsample4","IMsample5","IMsample6")
#log 2 of the DGE list 
log2_cpm_from_DGElist<-cpm(raw_DGEList, log=TRUE)
#Data frame from the log2 data frames include gene IDs and sample names
log_cpm_DF<-as_tibble(log2_cpm_from_DGElist, rownames= "geneID")
#get the colum names back 
colnames(log_cpm_DF)<-c("geneID",sample_names2)
log_cpm_DF
#Tidy data philosophy: each variable forms a column, each observation forms a row 
log_cpm_DF_pivot<-pivot_longer(log_cpm_DF, cols= AMsample1:IMsample6,names_to = "samples",values_to = "expression")
log_cpm_DF_pivot#for some reason the pivot_longer function do not work with the structure of the names AM-TB1 so I had change them to AMsample1,AMsample2,AMsample3...
#normalization of my data
DGElist_filtered_norm<-calcNormFactors(DGElist.filtered, method="TMM")
DGElist_filtered_norm
Log2_cpm_filtered_norm<-cpm(DGElist_filtered_norm, log=TRUE)
Log2_cpm_filtered_norm_df<-as_tibble(Log2_cpm_filtered_norm, rownames="geneID")
colnames(Log2_cpm_filtered_norm_df)<-c("geneID",sample_names2)
Log2_cpm_filtered_norm_df_pivot<-pivot_longer(Log2_cpm_filtered_norm_df,cols= AMsample1:IMsample6,names_to = "samples",values_to = "expression")
# use dplyr 'mutate' function to add new columns based on existing data
Log2_cpm_filtered_norm_df
mynewdata_df<-Log2_cpm_filtered_norm_df%>%
  mutate(AM_TB_negative = (AMsample1 + AMsample2 + AMsample3)/3,
         AM_TB_positive = (AMsample4 + AMsample5 + AMsample6 )/3,
         IM_TB_negative = (IMsample1 + IMsample2 + IMsample3)/3,
         IM_TB_positive = (IMsample4 + IMsample5 + IMsample6)/3,
         #now make columns comparing each of the averages above that you're interested in
         #LogFC means log of fold change 
         LogFC_AM_IM_neg = (AM_TB_negative - IM_TB_negative),
         LogFC_AM_IM_pos = (AM_TB_positive - IM_TB_positive),
         LogFC_AM_neg_pos = (AM_TB_positive - AM_TB_negative),
         LogFC_IM_neg_pos = (IM_TB_positive - IM_TB_negative)) %>% 
  mutate_if(is.numeric, round, 2)
# Differential Gene Expression -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 
#set up a design matrix 
group3<-factor(study_design$group)
deg<-model.matrix(~0 + group3)
colnames(deg)<-levels(group3)
# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v_DEGList_filtered_norm <- voom(DGElist_filtered_norm, deg, plot = TRUE)
#the stabilize variance of our model, apply precision weights to each gene 
#lowly expressed genes tend to have a higher stdv 
Fit<-lmFit(v_DEGList_filtered_norm, deg)
# Contrast matrix ----
contrast_matrix <- makeContrasts(infection_in_AM = AM_TB_Pos - AM_TB_Neg,
                                 infection_in_IM = IM_TB_Pos - IM_TB_Neg,
                                 Dexpression_celltypes=AM_TB_Neg - IM_TB_Neg, 
                                 Dexpression_celltypes_TB=AM_TB_Pos - IM_TB_Pos,
                                 levels=deg)


#Extract linear model fit 
Fits<-contrasts.fit(Fit,contrast_matrix )
#get bayesian stats for the linear model, this statistics will help identify the genes that are differentially express 
ebfit<- eBayes(Fits)

# TopTable to view DEGs -----
TopHits_AMs <- topTable(ebfit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
TopHits_IMs<- topTable(ebfit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
TopHits_AMs_IMs<- topTable(ebfit, adjust ="BH", coef=4, number=40000, sort.by="logFC")
#the toptable function from Limma select the top rank genes based on the contrast matrix 
#for each  pairwise comparison we are doing multiple statistical comparison tests 
#so we perform benjamin hershberger multiple testing corrections to reduce errors 
# convert to a tibble
TopHits_AMs_IMs_df <- TopHits_AMs_IMs %>%
  as_tibble(rownames = "geneID")
#negative fold changes indicates higher expression on IMs infected 
#positive fold changes indicates higher expression on AMs infected 

# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1 
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)
# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
vplot_AMs_IMs <- ggplot(TopHits_AMs_IMs_df ) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2)+
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot, TB infection in Alveolar Macrophages vs. TB infection in IMs",
       subtitle = "Up Regulation in AMs to the right, Up regulation in IMs to the left ",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot_AMs_IMs)
#on the x-axis we represent the magnitude of change in our pairwise comparison and in the y-axis 
#we use -log(adjpvalue) first to give the volcano shape , the genes on the top right and top left needed the least adjuste

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebfit, method="global", adjust.method="BH", p.value=0.01, lfc=1)

#adjust the adjusted p value and the lfc keep it on 1 
#essentially we are drawing a line right here and it is arbitrary there are some common 
#cutoffs such as a adjusted p.value=0.01 or p.value=0.05, and a logFoldChange=2 that is a 
#fold change of 4 


# take a look at what the results of decideTests looks like
#head(results)
#summary(results)
#vennDiagram(results, include="up")
#this is not useful with my comparison 

# retrieve expression data for your DEGs ----
head(v_DEGList_filtered_norm$E)
colnames(v_DEGList_filtered_norm$E) <- sample_id

differentGenes <- v_DEGList_filtered_norm$E[results[,4] !=0,]

#we keep all the genes that were other than 0 in the result matrix, all the genes that threw a 
#a significant difference in differential gene expression 
#check_genes

#convert your DEGs to a dataframe using as_tibble
differentGenes.df <- as_tibble(differentGenes, rownames = "geneID")
datatable(differentGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in AMs and IMs infected with TB',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:13), digits=2)

head(differentGenes)
dim(differentGenes)
# Creating Heatmaps  -----------
#this script creates heatmaps from your differentially expressed genes or transcripts
#and selects modules of co-expressed genes based on pearson correlations
library(tidyverse)
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3
# Choose your color pallette ----

# lots of easy ways to assemble your own color palette, including:
# 1). use 'colorRampPalette' function from the grDevices package
myheatcolors2 <- colorRampPalette(colors=c("purple","white","black"))(200)

# Cluster DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
# we use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(differentGenes), method="pearson")), method="complete") 
# hierarchical clustering is a type of unsupervised clustering. Related methods include K-means, SOM, etc 
# unsupervised methods are blind to sample/group identity
# in contrast, supervised methods 'train' on a set of labeled data.  
# supervised clustering methods include random forest, and artificial neural networks

#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(differentGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
#Cut the resulting tree and create color vector for clusters.  
#Vary the cut height (h =) to give more or fewer clusters, or use force k= number of clusters
#we'll look at these clusters in more detail later
module_assign <- cutree(clustRows, k=4)

#now assign a color to each module (makes it easy to identify and manipulate)
module_color <- rainbow(length(unique(module_assign)), start=0.1, end=0.9) 
module_color <- module_color[as.vector(module_assign)] 

heatmap.2(differentGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module_color,
          col=rev(myheatcolors2), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20),
          key.title = "1854 genes contrast between infected AMs and IMs") 

heatmaply(differentGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module_color,
          showticklabels=c(TRUE,FALSE),
          scale='row',
          key.title = "1854 genes contrast between infected AMs and IMs")


# OPTIONAL: simplify heatmap(not needed) ----
#notice that the heatmap includes ALL the columns from your dataset
#a useful way to simplify heatmaps, especially when there are many conditions, is to average your biological replicates and display only one column per condition
#rerun the heatmap script above using diffData.AVG as input instead of diffData
#colnames(differentGenes) <- study_design$group
colnames(differentGenes_celltypes) <- study_design$group

# Functional Enrichment Analysis  ----
# for the purposes of this script we'll want several data objects generated in previous scripts, including:
# 1) your normalized filtered expression data, in the form of a data matrix with symbols as rownames.
# 2) your study design file
# 3) your contrast matrix that lays out the pairwise comparisons you're interested in testing
# 4) Individual signatures or 'collections' of signatures to test for enrichment in your data.
# These signatures can be downloaded from gene signature databases such as MSigDB
# Signatures can also be custom made based on your interests.
# Signatures can also be pulled from R/Bioconductor as described below

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis

TopHits_AMs_IMs<- topTable(ebfit, adjust ="BH", coef=4, number=400, sort.by="logFC")
TopHits_AMs_IMs_subset <- TopHits_AMs_IMs[TopHits_AMs_IMs$logFC>=1,]
TopHits_IMs_AMs_subset <-TopHits_AMs_IMs[TopHits_AMs_IMs$logFC<=0,]
gost.res.AMs.up <- gost(rownames(TopHits_AMs_IMs_subset), organism = "mmusculus", correction_method = "fdr")
gost.res.IMs.up<- gost(rownames(TopHits_IMs_AMs_subset), organism = "mmusculus", correction_method = "fdr")
GOplotAMs_up<-gostplot(gost.res.AMs.up, interactive = F, capped = T) #set interactive=FALSE to get plot for publications
GOplotAM_up_interactive<-gostplot(gost.res.AMs.up, interactive = T, capped = T)

publish_gostplot(
  GOplotAMs_up, 
  highlight_terms = c("GO:0000070","GO:0007052","GO:00110564","GO:0022402","GO:0048285","GO:0051276","GO:0051726",
                      "GO:0098813","GO:0140014","GO:1902850"),
  filename = NULL,
  width = 10,
  height = 10)
#up regulated pathways in tb+ AMs 
GOplotIM_up_interactive<-gostplot(gost.res.IMs.up, interactive = T, capped = T)
#up regulated pathways in tb+ AMs

# option2: use the msigdb package to access up-to-date collections for Gene Set Enrichment Analysis----
# this option has the additional advantage of providing access to species-specific collections
# are also retrieved as tibbles
newdata_df<-Log2_cpm_filtered_norm_df%>%
  mutate(AM_TB_negative = (AMsample1 + AMsample2 + AMsample3)/3,
         AM_TB_positive = (AMsample4 + AMsample5 + AMsample6 )/3,
         IM_TB_negative = (IMsample1 + IMsample2 + IMsample3)/3,
         IM_TB_positive = (IMsample4 + IMsample5 + IMsample6)/3,
         #now make columns comparing each of the averages above that you're interested in
         #LogFC means log of fold change 
         LogFC_AM_IM_neg = (AM_TB_negative - IM_TB_negative),
         LogFC_AM_IM_pos = (AM_TB_positive - IM_TB_positive),
         LogFC_AM_neg_pos = (AM_TB_positive - AM_TB_negative),
         LogFC_IM_neg_pos = (IM_TB_positive - IM_TB_negative)) %>% 
  mutate_if(is.numeric, round, 2)
msigdbr_species()
hs_gsea <- msigdbr(species = "Mus musculus") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)%>%
  print(n=23)
# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
hs_gsea_h <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                     category = "H") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mynewdata_df_sub <- dplyr::select(newdata_df, geneID, LogFC_AM_IM_pos)
# construct a named vector
mydata_AMs_IM_gsea <- mynewdata_df_sub$LogFC_AM_IM_pos
names(mydata_AMs_IM_gsea) <- as.character(mynewdata_df_sub$geneID)
mydata_AMs_IMs_gsea <- sort(mydata_AMs_IM_gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA_res_AMs_IMs <- GSEA(mydata_AMs_IMs_gsea,pvalueCutoff = 0.05,pAdjustMethod = "BH", TERM2GENE=hs_gsea_h, verbose=FALSE)
myGSEA_res_AMs_IMs_df <- as_tibble(myGSEA_res_AMs_IMs@result)
#the @ symbol is used to assess a slot , it operates on f4 class objects in certain way simirlar to the $ operator 
myGSEA_AMs_IMs_table<-datatable(myGSEA_res_AMs_IMs_df, 
                            extensions = c('KeyTable', "FixedHeader"), 
                            caption = 'Signatures enriched in AMs infected with TB',
                            options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=10)
# create enrichment plots using the enrichplot package
plot_IMs_up1<-gseaplot2(myGSEA_res_AMs_IMs, 
                 geneSetID = c(3,4,19), #can choose multiple signatures to overlay in this plot
                 pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                 title = "upRegulated in IMs") #can also turn off this title
plot_AMs_up<-gseaplot2(myGSEA_res_AMs_IMs, 
                 geneSetID = c(11,15,9), #can choose multiple signatures to overlay in this plot
                 pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                 title = "upRegulated in AMs") #can also turn off this title
