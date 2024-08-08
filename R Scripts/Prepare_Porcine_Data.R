# Title: Cyclosporine A Accelerates Neurorecovery Transcriptional Trajectory in a Swine Model of Diffuse Traumatic Brain Injury
# Code created/curated by: Oluwagbemisola (Gbemi) Aderibigbe 
# Research Advisors: Levi B Wood, PhD; Susan S. Margulies, PhD
# Goal: Filtering, Normalization, and Outlier Removal 

# Install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("compare")
install.packages("gridGraphics")
install.packages("grid")
install.packages("tidyverse")
install.packages("factoextra")
install.packages("data.table")
install.packages("countToFPKM")
install.packages("writexl")
install.packages("viridis")
install.packages("flextable")
install.packages("ropls")
install.packages("DESeq2")
install.packages("biomaRt")
install.packages("org.Ss.eg.db")
install.packages("AnnotationDbi")
install.packages("vsn")


# Important
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # helps R know that folder that file is saved into is working directory
options(java.parameters = "-Xmx8g")  # gives Java more memory

# Load packages
library("gplots")
library("dendextend")
library("reshape")
library("matrixStats")
library("readxl")
library("compare")
library("grid")
library("gridGraphics")
library("DESeq2")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("boot")
library("apeglm")
library("ggrepel")
library("tidyverse")
library("grid")
library("gridExtra")
library("factoextra")
library("matrixStats")
library("ropls")
library("ggsci")
library("rio")
library ("readxl")
library("pacman")
library("biomaRt")
library("org.Ss.eg.db")
library("vsn")
library("countToFPKM")
library("tidyr")
library("data.table")
library("reshape2")
library("writexl")
library("viridis")
library("ggvenn")
library("flextable")
library("psych")
library("xlsx")

# Load source functions 
# Load WoodLabFunctions (https://github.com/afpybus/WoodLabFunctions/blob/main/R/rotate_opls.R)


### 1) CLEAN UP DATA 

# Import data
mRNA=read_excel("Sample Read Counts.xlsx")

# Specify characteristics of samples
Characteristics <- read_excel("Characteristics.xlsx") %>%   
  mutate(Region = case_when(
    Region == "Frontal" ~ "F",
    Region == "Hippocampus" ~ "Hip"
  )) %>%
  mutate(sampleID = str_c(Group,Number))

# Turn mRNA into row names
New_mRNA = as.matrix(mRNA[,2:ncol(mRNA)])
rownames(New_mRNA) = mRNA$id

# Filter data
iFilt=matrix(NA, ncol=1, nrow=nrow(New_mRNA))
for (i in 1:nrow(New_mRNA)) {
  if (sum(New_mRNA[i,]>=50)>20){iFilt[i]=TRUE}  # filter out genes with less than 50 counts in more than 20 samples
}
New_mRNA=New_mRNA[which(iFilt),]

# Change from ensembl ids to gene names
ensemblsIDS=as.data.frame(New_mRNA)
listEnsembl()
ensembl=useEnsembl(biomart="genes")
datasets=listDatasets(ensembl)
ensembl.con=useMart("ensembl", dataset="sscrofa_gene_ensembl")
attr=listAttributes(ensembl.con)
filters=listFilters(ensembl.con)

# Match ids
matched_ids=getBM(attributes= c('ensembl_gene_id','external_gene_name','start_position','end_position'),
                  filters= "ensembl_gene_id",
                  values=rownames(ensemblsIDS), 
                  mart=ensembl.con)

# Remove ensembl ids with no matched gene names 
matched_ids=matched_ids[!(is.na(matched_ids$external_gene_name) | matched_ids$external_gene_name==""), ]

# Merge corresponding ensembl id to gene names
setDT(ensemblsIDS, keep.rownames = "ensembl_gene_id")  #Make ensembl id a new column
NN=merge(ensemblsIDS, matched_ids[, c("ensembl_gene_id", "external_gene_name")], by="ensembl_gene_id")

# Make a new data frame with gene names as row names
NN2=NN[,-126]
NN3=as.matrix(NN2[,-1])
rownames(NN3)=NN$external_gene_name  # NN3 is your data with gene names 

# Group each sample id
MyData = tibble(.rows = ncol(NN3)) # create new tibble data frame
MyData$sampleID = colnames(NN3) # create a column with sample names in same order as data
MyData = left_join(MyData, Characteristics,by="sampleID") # add the experiment metadata
MyData = MyData %>%
  mutate(Condition = paste0(Type," ")) # create a condition variable with group
MyData$RegionandType= paste(MyData$Region, MyData$Type) #combine region and type



### 2) NORMALIZE DATA & PREPARE DATA FOR WGCNA 

# i) Normalize data with DESeq2
pheno1=MyData$RegionandType
condition = pheno1 
countData = NN3     
colData=(metaData<-data.frame(row.names = colnames(NN3), condition = pheno1))
colData$condition = factor(colData$condition) 
deseq.dat<- DESeqDataSetFromMatrix(countData,
                                   colData,
                                   design = ~ condition)

# Estimate size factor and normalize 
dds = estimateSizeFactors(deseq.dat)
sizeFactors(dds) # gives the size factor
New_mRNA_ = counts(dds, normalized=TRUE)  # gives DESeq2 normalized read counts


# ii) Prepare data for WGCNA using DESeq2 (WGCNA execution -- "WGCNA_Porcine_Analysis.R")
# Conduct variance stabilizing transformation
dds_s = DESeq(deseq.dat)
vsd = getVarianceStabilizedData(dds_s) 
datExpr = t(vsd)
export(datExpr,row.names=T, "datExpr.xlsx")



### 3) SORT & Z SCORE DATA FOR VISUALIZATION

# i) Select region and group
ideal_sorted = arrange(MyData,Region,Type)
sort_ind = match(ideal_sorted$sampleID,MyData$sampleID)
ROI = c("Hip", "F") # select region(s)
GOI = c("SHAM_4","Vehicle","Cyclosporine","Vehicle 7 day") # select group(s)
pdf = FALSE

# Sort by selected region and group
sortby <- sort_ind[sort_ind %in% which(
  MyData$Type %in% GOI &
    MyData$Region %in% ROI)]  
New_mRNA_sort <- New_mRNA_[,sortby]
geneID_analysis = mRNA$id
MyData2 <- MyData[sortby,]


# ii) Z score sorted data
FPKM_Z <- New_mRNA_sort %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(FPKM_Z) = colnames(New_mRNA_sort)
FPKM_Z$uniqueID = rownames(New_mRNA_sort)
FPKM_Z <- drop_na(FPKM_Z)   # remove rows with NA values
FPKM_Z_mat <- as.matrix(FPKM_Z[,1:{ncol(FPKM_Z)-1}])
rownames(FPKM_Z_mat) <- str_to_upper(FPKM_Z$uniqueID) # z scored data

# Re-arrange & Sort data (This is for making plots)
indSHAMF = which (MyData2$RegionandType == "F SHAM_4")
indInjured24F = which (MyData2$RegionandType == "F Vehicle")
indInjured7F = which (MyData2$RegionandType == "F Vehicle 7 day")
indCyclosporineF = which (MyData2$RegionandType == "F Cyclosporine")

indSHAMHip = which (MyData2$RegionandType == "Hip SHAM_4")
indInjured24Hip = which (MyData2$RegionandType == "Hip Vehicle")
indInjured7Hip = which (MyData2$RegionandType == "Hip Vehicle 7 day")
indCyclosporineHip = which (MyData2$RegionandType == "Hip Cyclosporine")

indSortFH=c(indSHAMF,indSHAMHip,indInjured24F,indInjured24Hip,
            indInjured7F, indInjured7Hip, indCyclosporineF, indCyclosporineHip) # can be modified depending on goal

MyData2=MyData2[indSortFH,]
New_mRNA_sort =New_mRNA_sort[,indSortFH]  # re-arranged count data 
FPKM_Z_mat=FPKM_Z_mat[, indSortFH]  # re-arranged z scored data (FPKM)



### 4) PLOT PCA, PLOT VIOLIN PLOTS & IDENTIFY OUTLIERS (Supplementary Fig. 2)

# Load WoodLabFunctions (https://github.com/afpybus/WoodLabFunctions/blob/main/R/rotate_opls.R)

# i) Conduct PCA and get PC1 and PC2 axis data
oplsOut <- opls(x = t(FPKM_Z_mat), predI=2)
rotateOut <- rotate_opls(oplsOut,degrees = 0)
dataPlot <- MyData2 %>%
  mutate(PC1 = rotateOut$T1) %>%
  mutate(PC2 = rotateOut$T2)  

# Plot PCA (Supplementary Fig. 2)
dataPlot$RegionandType=factor(dataPlot$RegionandType, levels =  c("F SHAM_4", "Hip SHAM_4", "F Vehicle", "Hip Vehicle",
                                                                  "F Vehicle 7 day", "Hip Vehicle 7 day", "F Cyclosporine", "Hip Cyclosporine"))
b=scores.plot(T1=dataPlot$PC1,T2=dataPlot$PC2,color = dataPlot$RegionandType, analysis="PCA")
b+ scale_color_manual(labels = c("F SHAM", "Hip SHAM", "F Injured 24 hours", "Hip Injured 24 hours",
                                 "F Injured 7 days", "Hip Injured 7 days", "F Cyclosporine", "Hip Cyclosporine"),
                      values = c("dodgerblue3", "cornflowerblue","firebrick1", "hotpink1", "goldenrod3", "tan2", "orchid4", "mediumpurple3")) + 
  labs(x = "Scores on PC1 (50.54%)", y = "Scores on PC2 (6.27%)")


# ii) Plot Violin Plots (Supplementary Fig. 2)
ggplot(dataPlot, aes(x=RegionandType, y=PC1, fill=RegionandType)) +
  geom_violin(trim=FALSE)+
  labs(
    title = "Comparison of PC1 Scores by Region and Group",
    x = "Region and Group",
    y = "PC1 Scores")+
  scale_fill_manual(labels = c("F SHAM", "Hip SHAM", "F Injured 24 hours", "Hip Injured 24 hours",
                               "F Injured 7 days", "Hip Injured 7 days", "F Cyclosporine", "Hip Cyclosporine"),
                    values = c("dodgerblue3", "cornflowerblue","firebrick1", "hotpink1", "goldenrod3", "tan2", "orchid4", "mediumpurple3"))+
  scale_x_discrete(labels = c(
    "F SHAM_4" = "F 
   SHAM",
    "F Vehicle" = "F 
   Inj 24 hrs",
    "F Vehicle 7 day" = "F 
   Inj 7 days",
    "F Cyclosporine" = "F 
   Cyclo",
    "Hip SHAM_4" = "Hip 
   SHAM",
    "Hip Vehicle" = "Hip 
   Inj 24 hrs",
    "Hip Vehicle 7 day" = "Hip 
   Inj 7 days",
    "Hip Cyclosporine" = "Hip 
   Cyclo"))+
  stat_compare_means(method="wilcox",
                     comparisons=list(c(1,2),c(3,4), c(5,6), c(7,8), c(1,3), c(2,4), c(1,5), c(2,6), c(1,7),c(2,8)),
                     p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)+
  geom_boxplot(width=.1, outlier.shape=NA, fill="white")+
  geom_point(color="black", size=1, position = position_jitter(w=0.05))+ theme_classic() 


# iii) Identify outliers
Groupss= dataPlot[c("PC1" , "PC2")] # get PC1 and PC2 data from above (i)
rownames(Groupss)=dataPlot$sampleID
Group.center= colMeans(Groupss)
Group.cov= cov(Groupss)
rad=qchisq(p=0.95, df=ncol(Groupss))
rad=sqrt(rad)
ellipse <- car::ellipse(center = Group.center , shape = Group.cov , radius = rad ,
                        segments = 150 , draw = FALSE)

ellipse <- as.data.frame(ellipse)
colnames(ellipse) <- colnames(Groupss)

figure <- ggplot(Groupss , aes(x = PC1 , y = PC2)) +
  geom_point(size = 2) +
  geom_polygon(data = ellipse , fill = "orange" , color = "orange" , alpha = 0.5)+
  geom_point(aes(Group.center[1] , Group.center[2]) , size = 5 , color = "blue") +
  geom_text( aes(label = row.names(Groupss)) , hjust = 1 , vjust = -1.5 ,size = 2.5 ) +
  ylab("PC2 scores") + xlab("PC1 scores")+ggtitle("PCA for SHAM Samples")

figure # plot to identify samples outside threshold (i.e outliers)
# Ensure to remove outlier samples then run 2,3 & 4 again then proceed to mRNA_Porcine_Analysis.R
