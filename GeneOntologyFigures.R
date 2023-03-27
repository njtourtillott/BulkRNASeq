library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)

cts.file <- "STAR_Gene_Counts.csv"
cts <- read.table(cts.file,sep=",",header=TRUE,stringsAsFactors=FALSE)
any(duplicated(rownames(cts$Gene_ID)))
rownames(cts) <- cts$Gene_ID
cts$Gene_ID <- NULL
universe.entrez <- na.omit(mapIds(org.Hs.eg.db, rownames(cts), 'ENTREZID', 'SYMBOL'))

# Read In Lists
compUp <- read.csv('sampleUp.csv')
up.genes <- compUp$X0
comparison<- 'sample'
compDown <-read.csv('sampleDown.csv')
down.genes <-compDown$X0


# Up GO
up.entrez <- na.omit(mapIds(org.Hs.eg.db, up.genes, 'ENTREZID', 'SYMBOL'))
up.go.test <- enrichGO(up.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
pdf(file=paste(comparison,"_up_go.pdf",sep=""),width=12,height=14)
print(dotplot(up.go.test,showCategory=40))
dev.off()
write.csv(data.frame(up.go.test),(file = paste(comparison,"_up_go.csv",sep="")))

# Down GO
down.entrez <- na.omit(mapIds(org.Hs.eg.db, down.genes, 'ENTREZID', 'SYMBOL'))
down.go.test <- enrichGO(down.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
pdf(file=paste(comparison,"_down_go.pdf",sep=""),width=12,height=14)
print(dotplot(down.go.test,showCategory=40))
dev.off()
write.csv(data.frame(down.go.test),(file = paste(comparison,"_down_go.csv",sep="")))


