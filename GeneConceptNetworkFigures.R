library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)



# Load In data
resOrdered.df <- read.delim("sample.txt",sep="\t")
comp <- 'sample'

# Filter Data & run GSEA
log2.fds <- resOrdered.df$log2FoldChange
names(log2.fds) <- rownames(resOrdered.df)
gene.list <- na.omit(log2.fds)
gene.list = sort(gene.list, decreasing = TRUE)

gmt <- read.gmt("hallmark.gmt")
pgmt <-read.gmt("c2.cp.v2022.1.Hs.symbols.gmt")
egmt <- GSEA(gene.list, TERM2GENE=gmt,minGSSize =6,
             pvalueCutoff = 0.05,pAdjustMethod = "fdr")

# Gene Concept Network Plot
p1 <- cnetplot(egmt, foldChange = log2.fds, showCategory = 5) + ggtitle(paste0(comp, '_geneConceptNetwork'))
pdf(paste0(comp,"_gcn.pdf"),width=20,height=25)
print(p1)
dev.off()

# #GSEA Heatmap
# h1 <- heatplot(egmt, foldChange = log2.fds, showCategory = 5)
# pdf( "heat.pdf",width=30,height=12)
# print(h1)
# dev.off()
# 
# # enrichment map
# pt <- enrichplot::pairwise_termsim(egmt)
# m1 <- enrichplot::emapplot(pt)
# pdf('map.pdf',width=20,height=25)
# print(m1)
# dev.off()

