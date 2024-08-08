# Title: Cyclosporine A Accelerates Neurorecovery Transcriptional Trajectory in a Swine Model of Diffuse Traumatic Brain Injury
# Code created/curated by: Oluwagbemisola (Gbemi) Aderibigbe 
# Research Advisors: Levi B Wood, PhD; Susan S. Margulies, PhD
# Goal: Differential Expression Analysis and Visualization

# Install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("heatmap3")
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
install.packages("ComplexHeatmap")
install.packages("EnhancedVolcano")
install.packages("UpSetR")

# Important
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # helps R know that folder that file is saved into is working directory
options(java.parameters = "-Xmx8g")  # gives Java more memory

# Load packages
library("gplots")
library("dendextend")
library("reshape")
library("matrixStats")
library("readxl")
library("heatmap3")
library("compare")
library("gridGraphics")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("boot")
library("DESeq2")
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
library("EnhancedVolcano")
library("ComplexHeatmap")
library("tidyr")
library("data.table")
library("reshape2")
library("writexl")
library("viridis")
library("ggvenn")
library("flextable")
library("psych")
library("xlsx")
library("UpSetR")

# Load source functions 
# Load WoodLabFunctions (https://github.com/afpybus/WoodLabFunctions/blob/main/R/rotate_opls.R)

###### Prepare data by running Prepare_Porcine_Data.R (needs to be ran before running these)


### 1) Re-arrange & Sort data (This is for making plots)
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



### 2) PLOT HEATMAPS FOR ENTIRE DATA (Fig. 1b)

# Make heatmap for entire data using Complex Heatmap 
fa=MyData2$RegionandType
fa_col=c("F SHAM_4"= 8, 
         "Hip SHAM_4"=8,
         "F Vehicle"= 9,
         "Hip Vehicle"= 9,
         "F Vehicle 7 day"= 7,
         "Hip Vehicle 7 day"=7)

dend2 = cluster_within_group(FPKM_Z_mat, fa) # conduct within group clustering
Heatmap(FPKM_Z_mat, cluster_columns = dend2, column_split = 6,
        row_title = "cluster_within_group",
        top_annotation = HeatmapAnnotation(foo = fa, col = list(foo = fa_col))) # visualize heatmap 



### 3) DIFFERENTIAL EXPRESSION ANALYSIS & PLOT VOLCANO PLOTS (Fig. 1c-f) & (Fig. 6a-d)

# i) Run stats for DESeq2
dds_stats= DESeq(deseq.dat)
head(results(dds_stats))
summary(results(dds_stats))
results(dds_stats,c("condition","F Vehicle", "F SHAM_4")) # two groups for comparison
contrast_kd <-  c("condition","F Vehicle", "F SHAM_4")    # two groups for comparison
res_table <- results(dds_stats, contrast=contrast_kd, alpha = 0.05)

# Identify genes with padj<0.05
table(res_table$padj<0.05)
Genename=rownames(res_table)[which(res_table$padj<0.05)]

# Export tables
Statitics_Raw_Data <- res_table@listData %>% as_tibble() %>%
  mutate(mRNA = res_table@rownames)
export(Statitics_Raw_Data,"Statitics_Raw_Data.csv") # CSV containing log2FC and p value is exported

# ii) Plot volcano plots (Fig. 1c-f)
keyvals <- ifelse(res_table$log2FoldChange < -0.5 & res_table$log2FoldChange > -0.99999  & res_table$pvalue < 0.05, '#8F8F8F',
                  ifelse(res_table$log2FoldChange < -1  & res_table$pvalue < 0.05, '#66C2A5',
                         ifelse(res_table$log2FoldChange > 0.5 & res_table$log2FoldChange < 0.99999 & res_table$pvalue < 0.05, '#8F8F8F',
                                ifelse(res_table$log2FoldChange > 1  & res_table$pvalue < 0.05, '#E78AC3',       
                                       '#8F8F8F'))))
keyvals[is.na(keyvals)] <- '#8F8F8F'
names(keyvals)[keyvals == '#8F8F8F'] <- 'up 0.5 < log2FC < 1'
names(keyvals)[keyvals == '#8F8F8F'] <- 'no change'
names(keyvals)[keyvals == '#8F8F8F'] <- 'down -0.5 < log2FC < 1'
names(keyvals)[keyvals == '#E78AC3'] <- 'up log2FC > 1'
names(keyvals)[keyvals == '#66C2A5'] <- 'down log2FC < -1'


EnhancedVolcano(res_table,
                lab = rownames(res_table),
                labSize=4,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Injured 1 day vs Sham - Frontal Cortex',
                selectLab = rownames(res_table)[which(names(keyvals) %in% c('up log2FC > 1', 'down log2FC < -1'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                colCustom = keyvals,
                colAlpha = 1,
                xlim = c(-3,3),
                ylim = c(0,12.5),
                pCutoff = 5e-2,
                FCcutoff=0.5,
                gridlines.major = TRUE,
                gridlines.minor = TRUE,
                legendPosition = "right",
                legendLabSize =10,
                legendIconSize = 5,
                border = 'full')



### 4) Sort OUT THE 6 TEMPORAL PATTERNS & PLOT PIECHARTS (Fig. 2c,d)
### (EARLY, LATE, PERSISTENT, INTENSIFIED, DELAYED & LATE)

# i) Import list of all genes meeting threshold of p value < 0.05 & log2FC > 1 OR < -1
Inj24=read_excel("Fron Inj 24 D.xlsx")
Inj7=read_excel("Fron Inj 7 D.xlsx")
Inj24vsInj7_Up = read_excel("Fron Up Inj 24 Down Inj 7.xlsx")
Inj24vsInj7_Down = read_excel("Fron Down Inj 24 Up Inj 7.xlsx")

# Template
Significant=intersect(Inj24,Inj7) # sig in both 1 day and 1 week
Only7days= setdiff(Inj7, Significant) # distinct to 1 week 
Only24hours=setdiff(Inj24,Significant) # distinct to 1 day

# Sort
Persistent1=setdiff(Significant,Inj24vsInj7_Up) # remove those with sig. difference between 1 day and 1 week
Persistent=setdiff(Persistent1,Inj24vsInj7_Down) # PERSISTENT

Intensified=intersect(Significant, Inj24vsInj7_Up) # INTENSIFIED

Delayed=intersect(Only7days,Inj24vsInj7_Up) # DELAYED

Transient=intersect(Only24hours,Inj24vsInj7_Down) # TRANSIENT

LateOther1=setdiff(Only7days, Inj24vsInj7_Down) # remove those with sig. difference between 1 day and 1 week
Late=setdiff(LateOther1, Inj24vsInj7_Up) # LATE

EarlyOther1=setdiff(Only24hours,Inj24vsInj7_Up) # remove those with sig. difference between 1 day and 1 week
Early=setdiff(EarlyOther1, Inj24vsInj7_Down) # EARLY 

# Save xlsx copy of temporal patterns
write_xlsx(Persistent,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Persistent\\Frontal Persistent.xlsx")
write_xlsx(Intensified,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Intensified\\Frontal Intensified.xlsx")
write_xlsx(Delayed,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Delayed\\Frontal Delayed.xlsx")
write_xlsx(Transient,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Transient\\Frontal Transient.xlsx")
write_xlsx(Late,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Late\\Frontal Late.xlsx")
write_xlsx(Early,"DESEQ2 Statistics\\Written Paper\\Grouping\\Downregulated\\Early\\Frontal Early.xlsx")


# ii) Plot Piecharts for Temporal patterns (Fig. 2c,d)
GPTED= read_excel(file.choose()) # import data with a table of total number of DEGs from above analysis

# Plot pie chart of temporal patterns
ggplot(GPTED, aes(x = 2, y = Value, fill = fct_inorder(Group))) +
  geom_col(width=1, size=1, color="white") +
  coord_polar("y", start = 0)+
  geom_text(aes(label = Value),
            position = position_stack(vjust = 0.5), size=6) +
  labs(x = NULL, y = NULL, fill = NULL, 
       title = "Hippocampus + Amygdala Downregulated DEGs Distribution") +
  theme_void()+ 
  scale_fill_manual(values=c("royalblue","deeppink","purple","yellow","yellowgreen","slategray")) 



### 5) CELL-TYPE SPECIFIC ANALYSIS & PLOT DOT PLOT (Fig. 2e,f)

# i) Identify cell-type specific DEGs in Late group (Late DEGs were identified in 4i)
A=read_excel("Astrocyte Genes.xlsx")
Astrocyte=distinct(A) # removes duplicates
colnames(Astrocyte)= c("mRNA")
BothUpA=intersect(Late,Astrocyte)

M=read_excel("Microglia Genes.xlsx")
Microglia=distinct(M) # removes duplicates
colnames(Microglia)= c("mRNA")
BothUpM=intersect(Late,Microglia)

N=read_excel("Neuronal Genes.xlsx")
Neuron=distinct(N) # removes duplicates
colnames(Neuron)= c("mRNA")
BothUpN=intersect(Late,Neuron)

O=read_excel("Oligondendrocyte Genes.xlsx")
Oligondendrocyte=distinct(O) # removes duplicates
colnames(Oligondendrocyte)= c("mRNA")
BothUpO=intersect(Late,Oligondendrocyte)

E=read_excel("Endothelial Genes.xlsx")
Endothelial=distinct(E) # removes duplicates
colnames(Endothelial)= c("mRNA")
BothUpE=intersect(Late,Endothelial)


# ii) Plot Dot Plot (Count) (Fig. 2e,f)
GeneCells=read_excel("Genes by Cells_Grouping_Up.xlsx") # import data with a table of total number of DEGs in each cell type from above analysis
GeneCells_=GeneCells[-6,]
GeneCells.long<-melt(GeneCells_) #melt as long list

# Plot dot plot for each cell type in all temporal patterns
ggplot(GeneCells.long, aes(y= fct_rev(variable), x= Cells)) + 
  geom_count(aes(size = value, color=Cells)) + 
  theme_classic()+
  scale_colour_manual(values=c("royalblue","deeppink","purple","yellow","yellowgreen","slategray"))+
  scale_size(range = c(0.5,20)) +
  labs(title = "Upregulated Gene Count per Cell", )+ xlab("Cells")+ ylab("Groups")



### 6) IDENTIFY TOP 10 DEGS & PLOT HEATMAPS (Fig. 3a,b) & (Fig. 4b) & (Fig. 5d-g)

# Set color bar
color.bar = function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) # create color bar based on z scored data
barColors = colorpanel(length(breakBarColors)-1, 'blue','white','red')

pdf("Colorbar.pdf", width=2,height=7,pointsize = 14) # export color bar to pdf
color.bar(barColors,-1.5,1.5,nticks=5)
dev.off()

# Set colors for group and region
pheno=MyData2$Condition
intcolorBar <- pheno
typeColorBar <- pheno

for (i in 1:length(unique(pheno))){
  intcolorBar[intcolorBar==unique(pheno)[i]] <- mypal[i]
}

regionpal <- brewer.pal(8,"Greens")  # palette for region and set color for region
regionColorBar <- case_when(
  MyData2$Region == "F" ~ regionpal[4],
  MyData2$Region == "Hip" ~ regionpal[7],
  MyData2$Region == "T" ~ regionpal[3],
)

cbPalette <- c("#A6CEE3","#FB9A99" ,"#E31A1C","#FFB347") # palette for group and set color for group
palette(cbPalette)
typeColorBar <- case_when(
  MyData2$Type == "Cyclosporine" ~ cbPalette[4],
  MyData2$Type == "Vehicle" ~ cbPalette[3],
  MyData2$Type == "Vehicle 7 day" ~ cbPalette[2],
  MyData2$Type == "SHAM_4" ~ cbPalette[1]
)

# i) Identify top 10 Upregulated and Downregulated Genes
T4U = read_excel("Hip Inj 24 U.xlsx") # 1 day upregulated list arranged based on log2FC values
T4U=T4U[1:10,]

T4D= read_excel("Hip Inj 24 D.xlsx") # 1 day downregulated list arranged based on log2FC values
T4D=T4D[1:10,]

S7U=read_excel("Hip Inj 7 U.xlsx") # 1 week upregulated list arranged based on log2FC values
S7U=S7U[1:10,]

S7D=read_excel("Hip Inj 7 D.xlsx") # 1 week downegulated list arranged based on log2FC values
S7D=S7D[1:10,]

# Identify the Transient DEGs and export normalized counts (Can also change DEGs to Top 10s)
DEGslist=FPKM_Z_mat[c('VAMP2','DRD2','CD79A','MADCAM1','SIX3'),] # Transient Downregulated Genes

# ii) Heatmap for Transient DEGs or Top 10 DEGs (Fig. 3a,b) & (Fig. 4b) & (Fig. 5d-g)
heatmap3(DEGslist,
         ColSideColors =cbind(Group=c(typeColorBar), Region=c(regionColorBar)),
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5),
         Rowv=NA, Colv=NA,
         labRow = rownames(DEGslist),mar=c(10,10),cexRow =1,cexCol = 0.65, main="Top 10 Hip 24 hours Up and Down",
         highlightCell=data.frame(rep(1:dim(DEGslist)[1],each=(dim(DEGslist)[2]+1)),
                                  rep(1:(dim(DEGslist)[2]+1),times=dim(DEGslist)[1]),'black',0.5),)



### 7) PLOT DOT PLOTS FOR ENRICHED GO TERMS (Fig. 3c) & (Fig. 4a,c) & (Fig. 5c) & (Fig. 6f-h)

# Import file from PANTHER (https://pantherdb.org/)
GOTerms_=read_excel(file.choose())
GOTerms=GOTerms_[-1:-6,c(-2,-3,-4,-5,-7)]  # remove all row and columns that are not needed
colnames(GOTerms)=c("Gobiologicalprocesscomplete","EnrichmentScore","pvalueFDR")

# Plot dot plot
ggplot(GOTerms, aes(y= reorder(Gobiologicalprocesscomplete, EnrichmentScore), x= log2(EnrichmentScore))) + 
  geom_count(aes(size = Count,color= pvalueFDR)) + 
  theme_bw()+
  scale_size(range = c(5,10))+
  labs(title = "GO Biological Processes", )+ xlab("Log2 Fold Enrichment")+ ylab("GO Biological Processes")+
  scale_colour_gradient(low = "blue",
                        high = "red")



### 8) PLOT BAR CHART FOR AXONAL INJURY (Fig. 4d)

# Import data
AxonalInj=read_excel("% Axonal Injury For Histogram Only Utilized Samples.xlsx")

# Plot Barplot
level_group=c("Sham", "Injured 24 hours", "Injured 7 days") # re-arrange data
ggbarplot(AxonalInj, x = "Group", y = "PercentAxonalInjury", 
          add = c("mean_se"),
          fill = "Group", palette = c("#A6CEE3","#E31A1C","#FB9A99"),
          position = position_dodge(0.9))+
  ggtitle("Percent Axonal Injury")+ylab("% Axonal Injury")+ xlab("Group")+
  geom_jitter(color = "black", width = 0.02)+
  stat_compare_means(method="wilcox",comparisons=list(c(1,2), c(2,3), c(1,3)),p.adjust.method = "bonferroni", label="p.format")



### 9) PLOT UPSET PLOTS (Fig. 5a,b)

# Get number of DEGS that overlap using "intersect" function
input <- c(
  Frontal.24 = 157,
  Hip.24 = 444,
  Frontal.7 = 520,
  Hip.7 = 1568,
  "Frontal.24&Hip.24"=127,
  "Frontal.7&Hip.7"=463,
  "Frontal.24&Frontal.7"=133,
  "Hip.24&Hip.7"=431,
  "Frontal.24&Hip.24&Frontal.7&Hip.7"=113)

# Plot
upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 4, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.9, 
      point.size = 2.8, 
      line.size = 1,
      mainbar.y.max = 3000
)



### 10) PLOT VENN DIAGRAMS (Fig. 6a-d)

# Import data 
All_Frontal4=read_excel("Fron Inj 24 D.xlsx")
All_Hip=read_excel("Fron Inj 7 D.xlsx")
All_Frontal3=read_excel("F Cyclosporine vs Sham Up in Sham.xlsx")

ALL_F4= as.vector(unlist(All_Frontal4))
ALL_H=as.vector(unlist(All_Hip))
ALL_F3=as.vector(unlist(All_Frontal3))
list_venn = list(Twentyfour=ALL_F4, Oneweek=ALL_H, Cyclosporine=ALL_F3)

# Plot Venn
ggvenn(list_venn, c("Twentyfour", "Oneweek", "Cyclosporine"),
       fill_color = c("#440154ff", "#21908dff","#D55E00"),
       stroke_size = 0.9, set_name_size = 5,
       text_size = 4, text_color = "black",  stroke_color = "black", stroke_linetype = "solid")
