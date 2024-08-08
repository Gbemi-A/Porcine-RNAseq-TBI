# Title: Cyclosporine A Accelerates Neurorecovery Transcriptional Trajectory in a Swine Model of Diffuse Traumatic Brain Injury
# Code created: Oluwagbemisola (Gbemi) Aderibigbe 
# Research Advisors: Levi B Wood, PhD; Susan S. Margulies, PhD
# Goal: WGCNA Analysis

#### IMPORTANT: Load Preliminary WGCNA codes from Alyssa Pybus (https://github.com/afpybus/bulk-RNAseq-rmTBI/blob/main/WGCNA.R)


### 1) Run WGCNA

#datExpr can be retrieved from Prepare_Porcine_Data.R

net = blockwiseModules(datExpr, power = 12, networkType = "signed",
                       deepSplit = 4,
                       minModuleSize = 10,
                       TOMDenom = "mean",
                       corType = "bicor",
                       mergeCutHeight = 0.15,
                       reassignThreshold = 0.05,
                       numericLabels = T, verbose = 3, maxBlockSize = 14000)


saveRDS(net,"DESEQ2 Statistics/Written Paper/WGCNA Output/net.RDS")



#### 2) PLOT THE DENDROGRAM (Fig. 7a)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    dendroLabels = FALSE,hang=0.03,addGuide = TRUE,guideHang = 0.05)



### 3) PLOT BOX PLOTS (Fig. 7b-d) 
## Ensure to sort data using (1) in mRNA_Porcine_Analysis.R
pheno=MyData2$Type
RawwME_data=MEs_sort %>%
  t() %>% 
  as_tibble() %>%
  mutate(Group =c("SHAM","SHAM","SHAM","SHAM","SHAM","SHAM","SHAM","SHAM",
                  "Injured 24 hours","Injured 24 hours","Injured 24 hours","Injured 24 hours","Injured 24 hours",
                  "Injured 24 hours","Injured 24 hours","Injured 24 hours","Injured 24 hours",
                  "Injured 7 days","Injured 7 days","Injured 7 days","Injured 7 days","Injured 7 days","Injured 7 days","Injured 7 days",
                  "Cyclosporine","Cyclosporine","Cyclosporine","Cyclosporine","Cyclosporine",
                  "Cyclosporine","Cyclosporine","Cyclosporine","Cyclosporine","Cyclosporine"
  ))%>%  
  mutate(Sample = colnames(MEs_sort))  

# Plot box plot for Module 10
pd = position_dodge(width = 0.5)
level_group=c("SHAM", "Injured 24 hours", "Cyclosporine", "Injured 7 days")
ggplot(RawwME_data, aes(x = factor(Group,level=level_group), y = `MEpurple`, fill = Group)) +
  stat_boxplot(geom="errorbar", position=pd, width=0.2) +
  geom_boxplot(width=0.5, position=pd)+ geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize=0.5, fill = "white")+
  ggtitle("Frontal Cortex ME 10 Score")+ylab("Read Count")+ xlab("Group")+
  stat_compare_means(method="wilcox",comparisons=list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4)),p.adjust.method = "bonferroni",label="p.format")+
  scale_fill_manual(labels = c("SHAM", "Injured 24 hours", "Cyclosporine", "Injured 7 days"),
                    breaks=c("SHAM", "Injured 24 hours", "Cyclosporine", "Injured 7 days"),
                    values = c("purple","purple","purple","purple"))+ theme_classic()



### 4) PLOT BAR PLOTS FOR CELLS (Supplementary Fig. 1) 

# Import file with number of cell type specific cells in each Module. Here is for Neurons
MEbarp= read_excel("NeuronME.xlsx")

barplot(height=MEbarp$Value, names=MEbarp$ME, 
        col=c("turquoise","blue","brown","yellow","green","red","black",
              "pink","magenta","purple","greenyellow","tan","salmon"),
        cex.names=0.7,
        xlab="Modules", 
        ylab="Number of genes", 
        main="Neurons",
        ylim=c(0,700)
)
