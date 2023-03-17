#Import the excell file from excell data sheet Filter and Normalize data----

library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
raw<-read.csv("lung_rawdata.csv",header = T)#import the data from the excell file 
genenames=raw[,1]#creating a list of gene names 
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
#sample_names<-study_design[,1]
sample_names2<-c("AMsample1","AMsample2","AMsample3","AMsample4","AMsample5","AMsample6","IMsample1","IMsample2","IMsample3","IMsample4","IMsample5","IMsample6")
#sample_names3<-c("AM_TB1","AM_TB2","AM_TB3","AM+TB1","AM+TB2","AM+TB3","IM_TB1","IM_TB2","IM_TB3","IM+TB1","IM+TB2","IM+TB3")
#log 2 of the DGE list 
log2_cpm_from_DGElist<-cpm(raw_DGEList, log=TRUE)
log2_cpm_from_DGElist
#Data frame from the log2 data frames include gene IDs and sample names
log_cpm_DF<-as_tibble(log2_cpm_from_DGElist, rownames= "geneID")
log_cpm_DF
#get the colum names back 
colnames(log_cpm_DF)<-c("geneID",sample_names2)
log_cpm_DF
#Tidy data philosophy: each variable forms a column, each observation forms a row 
log_cpm_DF_pivot<-pivot_longer(log_cpm_DF, cols= AMsample1:IMsample6,names_to = "samples",values_to = "expression")
log_cpm_DF_pivot#for some reason the pivot_longer function do not work with the structure of the names AM-TB1 so I had change them to AMsample1,AMsample2,AMsample3...

g1<-ggplot(log_cpm_DF_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
  
table(rowSums(raw_DGEList$counts==0)==12)
tokeep<-rowSums(cpm_from_rawDGElist>1)>=3#we filtered out all the counts that were equal to 0 
DGElist.filtered<-raw_DGEList[tokeep,]#from the DGElist we keep all the rows that had a count greater than 1 
DGElist.filtered
table(rowSums(DGElist.filtered$counts==0)==12)# to check how many rows with sums=0 were left 
#getting the CPM. and the log of the filtered data 
log_cpm_filtered<-cpm(DGElist.filtered, log=TRUE)
log_cpm_filtered_DF<-as_tibble(log_cpm_filtered, rownames="geneID")
colnames(log_cpm_filtered_DF)<-c("geneID",sample_names2)
log2_cpm_filtered_df_pivot<-pivot_longer(log_cpm_filtered_DF, cols= AMsample1:IMsample6,names_to = "samples",values_to = "expression")
#graphing the filtered data 
g2<-ggplot(log2_cpm_filtered_df_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
  
#normalization of my data
DGElist_filtered_norm<-calcNormFactors(DGElist.filtered, method="TMM")
DGElist_filtered_norm
Log2_cpm_filtered_norm<-cpm(DGElist_filtered_norm, log=TRUE)
Log2_cpm_filtered_norm_df<-as_tibble(Log2_cpm_filtered_norm, rownames="geneID")
colnames(Log2_cpm_filtered_norm_df)<-c("geneID",sample_names2)
Log2_cpm_filtered_norm_df_pivot<-pivot_longer(Log2_cpm_filtered_norm_df,cols= AMsample1:IMsample6,names_to = "samples",values_to = "expression")

g3<-ggplot(Log2_cpm_filtered_norm_df_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
plot_grid(g1, g2, g3, labels = c('A', 'B', 'C'), label_size = 12)
# Hierarchical clustering ---------------
library(tidyverse)
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # a layered 'grammar of tables' ggplot but for tables 
#Identify variable of interest in the study design file 
study_design
group1<-study_design$cell_type
group1<-factor(group1)
group2<-study_design$group
group2<-factor(group2)
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
mydistance <- dist(t(Log2_cpm_filtered_norm), method = "canberra") # we need to transpose the data first , because we are interested in the distance between genes and not the distance between samples 
mydistance
#other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
#when looking at the comparison of the distances between samples the highest distance verctors were between the AM-TB and the IM+TB 
myclusters <- hclust(mydistance, method = "complete") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(myclusters, labels=sample_names2)
#LEFT/RIGHT ORIENTATION HAS NO MEANING , the important is the vertical orientation 

#Principal Component Analysis PCA -------
pca_ana<-prcomp(t(Log2_cpm_filtered_norm), scale.=F, retx=T)
ls(pca_ana)
summary(pca_ana)
pca_ana$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca_ana$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca_ana) # A screeplot is a standard way to view eigenvalues for each PCA
pc_var<-pca_ana$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc_per<-round(pc_var/sum(pc_var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc_per
# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca_ana_df <- as_tibble(pca_ana$x)
ggplot(pca_ana_df) +
  aes(x=PC1, y=PC2, label=sample_names2, color = interaction(group1, group2)) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc_per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
#66% of variance accounted by PC1 which is cell type 
#18% of variance accounted by PC2 which is infected cell with TB and non infected cells 
cellt_plot<-ggplot(pca_ana_df) +
  aes(x=PC1, y=PC2, label=sample_names2, color = group1) +
  geom_point(size=4) +
  #geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc_per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
#visualize the principal components groups by colors in the graph 
infected_plot<-ggplot(pca_ana_df) +
  aes(x=PC1, y=PC2, label=sample_names2, color = group2) +
  geom_point(size=4) +
  #geom_label() +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc_per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

plot_grid(cellt_plot,infected_plot, labels = c('Cell type grouping ', 'TB infection grouping'), label_size = 12)

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca_ana_df <- pca_ana$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sample_names2,
             group = group2)

pca_pivot <- pivot_longer(pca_ana_df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca_pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
# Use dplyr 'verbs' to modify our dataframe ----
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
# the goal of this script is to identify differentially expressed genes (DEGs) and differential transcript usage (DTU)
# you should already know which pairwise comparisons are most important to you
# whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
# if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
# instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 and 4 to identify genes based on log fold-change
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
#the toptable function from Limma select the top rank genes based on the contrast matrix 
#for each  pairwise comparison we are doing multiple statistical comparison tests 
#so we perform benjamin hershberger multiple testing corrections to reduce errors 
# convert to a tibble
TopHits_AMs_df <- TopHits_AMs %>%
  as_tibble(rownames = "geneID")
TopHits_AMs_df#positive fold changes indicates that there are a higher expression on TB positive
#AMs 
TopHits_IMs_df<-TopHits_IMs%>%
  as_tibble(rownames = "geneID")
#a negative fold change indicates that fewer genes are expressed on TB+ IMs
#gt(TopHits_AMs_df)
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1 
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)
# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
vplot_AMs <- ggplot(TopHits_AMs_df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2)+
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot, TB infection in Alveolar Macrophages",
       subtitle = "TB infection in Alveolar Macrophages",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Now make the volcano plot above interactive with plotly
ggplotly(vplot_AMs)
#to the right(positive values are more expressed on AMs with TB infections )
vplot_IMs<-ggplot(TopHits_IMs_df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2, color="blue")+
  labs(title="Volcano plot, TB infection in Interstitial Macrophages",
     subtitle = "TB infection in nterstitial Macrophages",
     caption=paste0("produced on ", Sys.time())) +
  theme_bw()
ggplotly(vplot_IMs)
#on the x-axis we represent the magnitude of change in our pairwise comparison and in the y-axis 
#we use -log(adjpvalue) first to give the volcano shape , the genes on the top right and top left needed the least adjuste

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebfit, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results2<-decideTests(ebfit, method="global", adjust.method="BY", p.value=0.01, lfc=1)
#
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

differentGenes <- v_DEGList_filtered_norm$E[results[,1] !=0|results[,2] !=0,]
dim(differentGenes)
#
differentGenes2<-v_DEGList_filtered_norm$E[results2[,1] !=0|results2[,2] !=0,]
dim(differentGenes2)
#
differentGenes_celltypes <- v_DEGList_filtered_norm$E[results[,1] !=0|results[,2] !=0|results[,3] !=0,]
dim(differentGenes_celltypes)
#differentGenes_celltypes_TB<-v_DEGList_filtered_norm$E[results[,1] !=0|results[,2] !=0|results[,3] !=0|results[,4] !=0,]
#dim(differentGenes_celltypes_TB)
#we keep all the genes that were other than 0 in the result matrix, all the genes that threw a 
#a significant difference in differential gene expression 
#check_genes
head(differentGenes)
dim(differentGenes)
#convert your DEGs to a dataframe using as_tibble
differentGenes.df <- as_tibble(differentGenes, rownames = "geneID")
datatable(differentGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in AMs and IMs infected with TB',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:13), digits=2)

#data frame including differentially expressed genes with a lfc=1 and adjust Method: BY
differentGenes.df2 <- as_tibble(differentGenes2, rownames = "geneID")
datatable(differentGenes.df2, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in AMs and IMs infected with lfc=2',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:13), digits=2)

# Creating Heatmaps  -----------
#this script creates heatmaps from your differentially expressed genes or transcripts
#and selects modules of co-expressed genes based on pearson correlations
library(tidyverse)
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(gameofthrones) #because...why not.  Install using 'devtools::install_github("aljrico/gameofthrones")'
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3
# Choose your color pallette ----
#Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcolors1 <- bluered(75) #this is from the 'colorpanel' function in gplots (same package that heatmap.2 comes from)
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

# lots of easy ways to assemble your own color palette, including:
# 1). use 'colorRampPalette' function from the grDevices package
myheatcolors2 <- colorRampPalette(colors=c("blue","white","brown"))(100)
# 2). use rcolorbrewer to choose any palette by name and n colors from that palette
myheatcolors3 <- brewer.pal(name="RdBu", n=11)
# 3). paste in your own hex codes using the Sip app (or other tools)
myheatcolors3 <- c("#fed976", "#268f9c")
#differentGenes_all <- v_DEGList_filtered_norm$E[results[,1] !=0|results[,2] !=0,|results[,2] !=0,]

# Cluster DEGs ----
#begin by clustering the genes (rows) in each set of differentially expressed genes
# we use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(differentGenes), method="pearson")), method="complete") 
clustRows_lfc1<-hclust(as.dist(1-cor(t(differentGenes2), method="pearson")), method="complete") 
# hierarchical clustering is a type of unsupervised clustering. Related methods include K-means, SOM, etc 
# unsupervised methods are blind to sample/group identity
# in contrast, supervised methods 'train' on a set of labeled data.  
# supervised clustering methods include random forest, and artificial neural networks

#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(differentGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
clustColumns_lfc1<-hclust(as.dist(1-cor(differentGenes2, method="spearman")), method="complete")
#note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes
#see lecture slides for an example using a mock dataset.

#Cut the resulting tree and create color vector for clusters.  
#Vary the cut height (h =) to give more or fewer clusters, or use force k= number of clusters
#we'll look at these clusters in more detail later
module_assign <- cutree(clustRows, k=4)
module_assign_lfc1<-cutree(clustRows_lfc1, k=4)

#now assign a color to each module (makes it easy to identify and manipulate)
module_color <- rainbow(length(unique(module_assign)), start=0.1, end=0.9) 
module_color <- module_color[as.vector(module_assign)] 
module_color_lfc1<-rainbow(length(unique(module_assign_lfc1)), start=0.1, end=0.9)
module_color_lfc1 <- module_color_lfc1[as.vector(module_assign_lfc1)]
heatmap.2(differentGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module_color,
          col=rev(myheatcolors2), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20),
          key.title = "1559genes_lfc=1") 

heatmaply(differentGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module_color,
          showticklabels=c(TRUE,FALSE),
          scale='row',
          key.title = "1559genes_lfc=1")
#heat map with all the 1164 genes same cluster comparison , adjustment method: BY  
heatmap.2(differentGenes2, 
          Rowv=as.dendrogram(clustRows_lfc1), 
          Colv=as.dendrogram(clustColumns_lfc1),
          RowSideColors=module_color_lfc1,
          col=rev(myheatcolors2), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20),
          key.title = "1164, AdM=BY") 
heatmaply(differentGenes2,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows_lfc1),
          RowSideColors=module_color_lfc1,
          showticklabels=c(TRUE,FALSE),
          scale='row',
          key.title = "1164,adM=BY")

# OPTIONAL: simplify heatmap(not needed) ----
#notice that the heatmap includes ALL the columns from your dataset
#a useful way to simplify heatmaps, especially when there are many conditions, is to average your biological replicates and display only one column per condition
#rerun the heatmap script above using diffData.AVG as input instead of diffData
#colnames(differentGenes) <- study_design$group
colnames(differentGenes_celltypes) <- study_design$group

#now an old function from the limma package to average your replicates 
#differentGenes.AVG <- avearrays(differentGenes)
#differrentGenes_celltypes.AVG<-avearrays(differentGenes_celltypes)

##alternatively, decide exactly which columns you want to show, and modify the heatmap accordingly
#this is how it would look using base R
#diffGenes.subset <- diffGenes[,c(1,4,7)]
##now repeat heatmap using only these selected columns
#clustRows2 <- hclust(as.dist(1-cor(t(differentGenes.AVG), method="pearson")), method="complete") 
#clustColumns2 <- hclust(as.dist(1-cor(differentGenes.AVG, method="spearman")), method="complete") 
#heatmap.2(differentGenes.AVG, 
          #Rowv=as.dendrogram(clustRows2), 
          #Colv=as.dendrogram(clustColumns2),
          #RowSideColors=module_color,
          #col=rev(myheatcolors2), scale='row', labRow=NA,
          #density.info="none", trace="none",  
          #cexRow=1, cexCol=1, margins=c(8,20)) 
#heatmaply(differentGenes.AVG,
          #colors = myheatcolors2,
          #Rowv=as.dendrogram(clustRows2),
          #RowSideColors=module_color,
          #showticklabels=c(TRUE,FALSE),
          #scale='row')

#simplify heat map for the 1444 genes 
#clustRows2_celltypes  <- hclust(as.dist(1-cor(t(differrentGenes_celltypes.AVG), method="pearson")), method="complete") 
#clustColumns2_celltypes <- hclust(as.dist(1-cor(differrentGenes_celltypes.AVG, method="spearman")), method="complete") 
#heatmap.2(differrentGenes_celltypes.AVG, 
          #Rowv=as.dendrogram(clustRows2_celltypes), 
          #Colv=as.dendrogram(clustColumns2_celltypes),
          #RowSideColors=module_color_celltypes,
          #col=rev(myheatcolors2), scale='row', labRow=NA,
          #density.info="none", trace="none",  
          #cexRow=1, cexCol=1, margins=c(8,20)) 
#heatmaply(differrentGenes_celltypes.AVG,
          #colors = myheatcolors2,
          #Rowv=as.dendrogram(clustRows),
          #RowSideColors=module_color,
          #showticklabels=c(TRUE,FALSE),
          #scale='row')

# View modules of co-regulated genes ----
# view your color assignments for the different clusters
names(module_color) <- names(module_assign) 

module_assign_df <- as_tibble(as.list(module_assign))
module_assign_pivot <- pivot_longer(module_assign_df, # dataframe to be pivoted
                                    cols = 1:1559, # column names to be stored as a SINGLE variable
                                    names_to = "geneID", # name of that new variable (column)
                                    values_to = "module") # name of new variable (column) storing all the values (data)

module_assign_pivot <- module_assign_pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "AM-",
    module == 2 ~ "IM+",
    module == 3 ~ "AM+",
    module == 4 ~ "IM-"))


ggplot(module_assign_pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
modulePick <-4 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)
#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- differentGenes[names(module_assign[module_assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=TRUE, 
          labRow = NA,
          col=rev(myheatcolors2), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module_color[module_assign%in%modulePick], margins=c(8,20)) 
dim(myModule)
#Module 2: 
# Include the genes that are  up regulated in infected IMs 398 genes of interest , 
#Module 1: 
#include genes that are downregulated in infected IMs 454 genes of interest 
#Module 3:
#genes upregulated in infected AMs 196 genes of interest  
#module4: 
#downregulated genes in infected, 511 genes of interest AMs 

# Export modules for downstream analysis ----
#prints out genes in the order you see them in the cluster
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- differentGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"downregulated genes infected AMs ")

modulePick_upAMs <-3
myModule_upAMs <- differentGenes[names(module_assign[module_assign %in% modulePick_upAMs]),]
modulePick_upIMs<-2
myModule_upIMs <- differentGenes[names(module_assign[module_assign %in% modulePick_upIMs]),]
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
TopHits_AMs <- topTable(ebfit, adjust ="BH", coef=1, number=250, sort.by="logFC")
TopHits_IMs<- topTable(ebfit, adjust ="BH", coef=2, number=250, sort.by="logFC")
TopHitsIMs_subset <- TopHits_IMs[TopHits_IMs$logFC>1,]
#coef=1 is referring to contrast matrix of comparisons between infected and non infected groups
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
TopHitsAMs_treshold<-rownames(TopHits_AMs[TopHits_AMs$logFC > 1,]) 
TopHitsAM_subset <- TopHits_AMs[TopHits_AMs$logFC>1,]
gost.res.AMs <- gost(rownames(TopHitsAM_subset), organism = "mmusculus", correction_method = "fdr")
gost.res.IMs <- gost(rownames(TopHitsIMs_subset), organism = "mmusculus", correction_method = "fdr")
gost.res.Module_upAMs<- gost(rownames(myModule_upAMs), organism = "mmusculus", correction_method = "fdr")
gost.res.Module_upIMs<- gost(rownames(myModule_upIMs), organism = "mmusculus", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
GOplot1<-gostplot(gost.res.AMs, interactive = F, capped = T) #set interactive=FALSE to get plot for publications
GOplotAM_interactive<-gostplot(gost.res.AMs, interactive = T, capped = T)
#GOplot1<-
GOplotAM_module<-gostplot(gost.res.Module_upAMs, interactive = T, capped = T)
# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
GoplotIM_topHits<-gostplot(gost.res.IMs, interactive = F, capped = T)
gostplot(gost.res.IMs, interactive = T, capped = T)
gostplot(gost.res.Module_upIMs, interactive = T, capped = T)
publish_gostplot(
  GOplot1, 
 highlight_terms = c("GO:0006954","GO:0051239","GO:0002682","GO:000695","GO:0002376","GO:0009605","GO:0051240",
 "GO:0071944","GO:0005615","GO:0005576","GO:0009986","GO:0005102","GO:0008289","GO:00035325","GO:0002682","KEGG:04060",
 "REAC:R-MMU-168256","CORUM:6922","CORUM:146","CORUM:6938","WP:WP4466","WP:WP4335"),
  filename = NULL,
  width = 10,
  height = 10)
#IMs enriched pathways  
publish_gostplot(
  GoplotIM_topHits, 
  highlight_terms = c("GO:0002376","GO:0006954","GO:0006952","GO:0006950","GO:0009605","GO:051239","GO:0005615","REAC:R-MMU-6798695",
  "REAC:R-MMU-164249","REAC:R-MMU-5621480", "MIRNA:mmu-miR-3100-3p"),
  filename = NULL,
  width = 10,
  height = 10)

# option2: use the msigdb package to access up-to-date collections for Gene Set Enrichment Analysis----
# this option has the additional advantage of providing access to species-specific collections
# are also retrieved as tibbles
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
mydata_df_sub <- dplyr::select(mynewdata_df, geneID, LogFC_AM_neg_pos,LogFC_IM_neg_pos)
# construct a named vector
mydata_AMs_gsea <- mydata_df_sub$LogFC_AM_neg_pos
names(mydata_AMs_gsea) <- as.character(mydata_df_sub$geneID)
mydata_AMs_gsea <- sort(mydata_AMs_gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA_res_AMs <- GSEA(mydata_AMs_gsea,pvalueCutoff = 0.05,pAdjustMethod = "BH", TERM2GENE=hs_gsea_h, verbose=FALSE)
myGSEA_res_AMs_df <- as_tibble(myGSEA_res_AMs@result)
#the @ symbol is used to assess a slot , it operates on f4 class objects in certain way simirlar to the $ operator 
myGSEA_AMs_table<-datatable(myGSEA_res_AMs_df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in AMs infected with TB',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=10)
# create enrichment plots using the enrichplot package
plot1<-gseaplot2(myGSEA_res_AMs, 
          geneSetID = c(2,4,7), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA_res_AMs$Description[47]) #can also turn off this title
plot2<-gseaplot2(myGSEA_res_AMs, 
          geneSetID = c(12,13,23), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA_res_AMs$Description[47]) #can also turn off this title
#Pulling out IMs information 
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
# construct a named vector
mydata_IMs_gsea <- mydata_df_sub$LogFC_IM_neg_pos
names(mydata_IMs_gsea) <- as.character(mydata_df_sub$geneID)
mydata_IMs_gsea <- sort(mydata_IMs_gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA_res_IMs <- GSEA(mydata_IMs_gsea, TERM2GENE=hs_gsea_h, verbose=FALSE)
myGSEA_res_IMs_df <- as_tibble(myGSEA_res_IMs@result)
#the @ symbol is used to assess a slot , it operates on f4 class objects in certain way simirlar to the $ operator 
myGSEA_IMs_table<-datatable(myGSEA_res_IMs_df, 
                            extensions = c('KeyTable', "FixedHeader"), 
                            caption = 'Signatures enriched in IMs infected with TB',
                            options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=10)
# create enrichment plots using the enrichplot package
plot1_IMs<-gseaplot2(myGSEA_res_IMs, 
                 geneSetID = c(3,4,5,6), #can choose multiple signatures to overlay in this plot
                 pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                 title = myGSEA_res_IMs$Description[47]) #can also turn off this title
plot2_IMs<-gseaplot2(myGSEA_res_IMs, 
                 geneSetID = c(9,11,14,16), #can choose multiple signatures to overlay in this plot
                 pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                 title = myGSEA_res_AMs$Description[47]) #can also turn off this title
# add a variable to this result that matches enrichment direction with phenotype
myGSEA_res_AMs_df <- myGSEA_res_AMs_df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))
#IMs results 
myGSEA_res_IMs_df <- myGSEA_res_IMs_df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))
# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA_res_AMs_df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="Enriched  pathways in Hallmark signature from MsigDB in AMs ")
  theme_bw()
#IMs results 
ggplot(myGSEA_res_IMs_df[1:26,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="gray", high="green") +
  labs(title="Enriched  pathways in Hallmark signature from MsigDB in IMs ")
  theme_bw()
#New ImmunoMSigDB----
hs_gsea_C7 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                     category = "C7") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol)
# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA_res_AMs_Immune <- GSEA(mydata_AMs_gsea,pvalueCutoff = 0.05,pAdjustMethod = "BH", TERM2GENE=hs_gsea_C7, verbose=FALSE)
myGSEA_res_AMs_Immune_df <- as_tibble(myGSEA_res_AMs_Immune@result)
myGSEA_AMs_table_Immune<-datatable(myGSEA_res_AMs_Immune_df, 
                            extensions = c('KeyTable', "FixedHeader"), 
                            caption = 'Signatures enriched in AMs infected with TB, ImmunoMSigDB',
                            options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
                            formatRound(columns=c(2:10), digits=10)
#Bubble plot 
myGSEA_res_AMs_Immune_df<- myGSEA_res_AMs_Immune_df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))
ggplot(myGSEA_res_AMs_Immune_df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="grey", high="purple") +
  labs(title="Enriched  pathways in Immunology C7 signature from MsigDB in AMs ")
  theme_bw()
#over represented pathways in Immunology C7 signature of MsigDB in AMs 
#getting IMs data 
myGSEA_res_IMs_Immune <- GSEA(mydata_IMs_gsea,pvalueCutoff = 0.05,pAdjustMethod = "BH", TERM2GENE=hs_gsea_C7, verbose=FALSE)
myGSEA_res_IMs_Immune_df <- as_tibble(myGSEA_res_IMs_Immune@result)
myGSEA_IMs_table_Immune<-datatable(myGSEA_res_IMs_Immune_df, 
                            extensions = c('KeyTable', "FixedHeader"), 
                            caption = 'Pathways enriched in AMs infected with TB, ImmunoMsigDB',
                            options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=10)
#Bubble plot 
myGSEA_res_IMs_Immune_df<- myGSEA_res_IMs_Immune_df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))
ggplot(myGSEA_res_IMs_Immune_df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="black", high="blue") +
  labs(title="Enriched pathways in Immunology C7 signature from MsigDB in IMs ")
  theme_bw()
# Competitive GSEA using CAMERA----
# for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially expressed as genes outside the set
# first let's create a few signatures to test in our enrichment analysis
mySig_fromTopHits <- rownames(TopHits_AMs) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig_random <- sample((rownames(v_DEGList_filtered_norm$E)), size = 250, replace = FALSE)
collection <- list(real = mySig_fromTopHits, fake = mySig_random)
# now test for enrichment using CAMERA
#Testing enrichment in the Top DE genes-list corresponding to AMs 
camera.result <- camera(v_DEGList_filtered_norm$E, collection, deg, contrast_matrix[,1]) 
camera_df <- as_tibble(camera.result, rownames = "setName")
camera_df
#Testing enrichment in the Top DE genes-list corresponding to IMs 
mySig_fromTopHits_IMs <- rownames(TopHits_IMs) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig_random_IMs <- sample((rownames(v_DEGList_filtered_norm$E)), size = 250, replace = FALSE)
collection <- list(real = mySig_fromTopHits_IMs, fake = mySig_random_IMs)
# now test for enrichment using CAMERA
camera.result_IMs <- camera(v_DEGList_filtered_norm$E, collection, deg, contrast_matrix[,2]) 
camera_df_IMs <- as_tibble(camera.result_IMs, rownames = "setName")
camera_df_IMs
#my top differentially expressed genes in IMs are down regulated 

