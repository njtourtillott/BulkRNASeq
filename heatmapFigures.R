library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)

## Load in Log2Fold Chnage data 
resOrdered.df <- read.delim("sample.txt",sep="\t")
samp <- 'sample' #enter Sample name here

## Load in Pathways of interest
# from gmt file
gmt <- read.gmt("/hallmark.gmt") 
pathway1 <- gmt[gmt$term == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]
# or from String-db
integrin <- read.delim('String_integrinMediatedSignalingPathway.tsv', sep = '\t')


## Create dataframe to plot
# needs structure of : Genes|Pathway|Log2FoldChange

# Choose genes of interest
genes <- c('ITGB3') # your genes here or read in if needed
genes <- data.frame(genes)
genes2 <- genes # for more than one pathway

genes$Pathway <- 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'
genes$log2FoldChange <- NA
# Check if genes in pathways of interest, if so insert LFC if not 0/NA
for(i in length(genes$genes)){
  if(is.element(genes$genes[i], pathway1$gene)){
   genes$log2FoldChange[i] <- resOrdered.df[genes$genes[i], 'log2FoldChange']
  }
} 

genes2$Pathway <- 'Integrin'
genes2$log2FoldChange <- NA
for(i in length(genes$genes)){
  if(is.element(genes$genes[i], integrin$X.node)){
    genes2$log2FoldChange[i] <- resOrdered.df[genes$genes[i], 'log2FoldChange']
  }
} 

# Create dataframe
df <- rbind(genes, genes2) # add extra pathways as needed

## Plot
hMap <- ggplot(df, aes(x = genes,
                          y = Pathway,
                          fill = log2FoldChange))+geom_tile()
pdf(paste0(samp, "_pathwayHeatmap.pdf"),width=20,height=15)
print(hMap)
dev.off()


