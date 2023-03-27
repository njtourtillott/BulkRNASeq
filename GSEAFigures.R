library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Read in GMT and Gene LFC data
gmt <- read.gmt("hallmark.gmt")
resOrdered.df <- read.delim("sample.txt",sep="\t")
comparison <- 'sample'

#read in overlap genes to list
compUp <- read.csv('UpRegulatedOverlap.csv')
up.genes <- compUp$X
compDown <-read.csv('10DayDownRegulatedOverlap.csv')
down.genes <-compDown$X
#combine up and down overlaps 
gene.overlap <- c(up.genes, down.genes)

  
# Filter Comparison.txt to just LFC of overlap genes
resOrdered.df <- resOrdered.df[rownames(resOrdered.df) %in% gene.overlap,]
write.table(resOrdered.df,file=paste(comparison,".txt",sep=""),sep="\t")
log2.fds <- resOrdered.df$log2FoldChange
names(log2.fds) <- rownames(resOrdered.df)
gene.list <- na.omit(log2.fds)
gene.list = sort(gene.list, decreasing = TRUE)


# Run GSEA on Gene list
egmt <- GSEA(gene.list, TERM2GENE=gmt,minGSSize =6,
             pvalueCutoff = 0.05,pAdjustMethod = "fdr")

# Plot GSEA #
pdf(paste0(comparison, "_gseaHallmark.pdf"),width=12,height=10)
print(ggplot(egmt@result, aes(reorder(Description, NES), NES)) +
        geom_col(aes(fill=p.adjust)) +
        coord_flip() +
        geom_text(aes(label = signif(p.adjust,digits=2)), size=5, nudge_y  = -0.2, colour = "black") +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title=comparison) +
        theme_minimal()+ theme(legend.text = element_text(colour="black", size=10),
                               axis.text.y =element_text(colour="black", size=15),
                               axis.title.y =element_text(colour="black", size=15),
                               axis.text.x =element_text(colour="black", size=15),
                               axis.title.x =element_text(colour="black", size=15))+
        scale_fill_gradient2(low = "red", high = "blue",mid = "purple",limits=c(0,1),midpoint = 0.5)
)
dev.off()
