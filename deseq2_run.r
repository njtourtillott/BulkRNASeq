library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)


comparisons <- list()
## Start
# 1
comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]] <- list()
comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]][["columns"]] <- c("PC9.Osi_24_1", "PC9.Osi_24_2", "PC9.OsiplusJH_24_1", 'PC9.OsiplusJH_24_2')

comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]][["condition"]] <- c("PC9.Osi_24", 'PC9.Osi_24', 'PC9.OsiplusJH_24','PC9.OsiplusJH_24' 
)

comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]][["batch"]] <- c("1","2",
                                                                              '1', '2')

comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]][["title"]] <- "Change in Expression\nPC9.OsiplusJH_24_vs_PC9.Osi_24"
comparisons[["PC9.OsiplusJH_24_vs_PC9.Osi_24"]][["contrast"]] <- c("condition","PC9.OsiplusJH_24", "PC9.Osi_24")

#2
comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]] <- list() 
comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]][["columns"]] <- c("X4006.Osi_24_1", "X4006.Osi_24_2", "X4006.OsiplusJH_24_1", 'X4006.OsiplusJH_24_2')

comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]][["condition"]] <- c("X4006.Osi_24", 'X4006.Osi_24', 'X4006.OsiplusJH_24','X4006.OsiplusJH_24' 
)

comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]][["batch"]] <- c("1","2",
                                                                '1', '2')

comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]][["title"]] <- "Change in Expression\nX4006.OsiplusJH_24_vs_X4006.Osi_24"
comparisons[["X4006.OsiplusJH_24_vs_X4006.Osi_24"]][["contrast"]] <- c("condition","X4006.OsiplusJH_24", "X4006.Osi_24")

#3
comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]] <- list()
comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]][["columns"]] <- c("PC9.Osi_72_1", "PC9.Osi_72_2", "PC9.OsiplusJH_72_1", 'PC9.OsiplusJH_72_2')

comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]][["condition"]] <- c("PC9.Osi_72", 'PC9.Osi_72', 'PC9.OsiplusJH_72','PC9.OsiplusJH_72' 
)

comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]][["batch"]] <- c("1","2",
                                                                '1', '2')

comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]][["title"]] <- "Change in Expression\nPC9.OsiplusJH_72_vs_PC9.Osi_72"
comparisons[["PC9.OsiplusJH_72_vs_PC9.Osi_72"]][["contrast"]] <- c("condition","PC9.OsiplusJH_72", "PC9.Osi_72")

#4
comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]] <- list() 
comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]][["columns"]] <- c("X4006.Osi_72_1", "X4006.Osi_72_2", "X4006.OsiplusJH_72_1", 'X4006.OsiplusJH_72_2')

comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]][["condition"]] <- c("X4006.Osi_72", 'X4006.Osi_72', 'X4006.OsiplusJH_72','X4006.OsiplusJH_72' 
)

comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]][["batch"]] <- c("1","2",
                                                                    '1', '2')

comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]][["title"]] <- "Change in Expression\nX4006.OsiplusJH_72_vs_X4006.Osi_72"
comparisons[["X4006.OsiplusJH_72_vs_X4006.Osi_72"]][["contrast"]] <- c("condition","X4006.OsiplusJH_72", "X4006.Osi_72")


#5
comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]] <- list()
comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]][["columns"]] <- c("PC9.Osi_10_day_or_DTC_1", "PC9.Osi_10_day_or_DTC_2", "PC9.Osi_plusJH_10_day_or_DTC_1", 'PC9.Osi_plusJH_10_day_or_DTC_2')

comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]][["condition"]] <- c("PC9.Osi_10_day_or_DTC", 'PC9.Osi_10_day_or_DTC', 'PC9.Osi_plusJH_10_day_or_DTC','PC9.Osi_plusJH_10_day_or_DTC' 
)

comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]][["batch"]] <- c("1","2",
                                                                '1', '2')

comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]][["title"]] <- "Change in Expression\nPC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"
comparisons[["PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC"]][["contrast"]] <- c("condition","PC9.Osi_plusJH_10_day_or_DTC", "PC9.Osi_10_day_or_DTC")

#6
comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]] <- list() 
comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]][["columns"]] <- c("X4006.Osi_10_day_or_DTC_1", "X4006.Osi_10_day_or_DTC_2", "X4006.Osi_plusJH_10_day_or_DTC_1", 'X4006.Osi_plusJH_10_day_or_DTC_2')

comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]][["condition"]] <- c("X4006.Osi_10_day_or_DTC", 'X4006.Osi_10_day_or_DTC', 'X4006.Osi_plusJH_10_day_or_DTC','X4006.Osi_plusJH_10_day_or_DTC' 
)

comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]][["batch"]] <- c("1","2",
                                                                    '1', '2')

comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]][["title"]] <- "Change in Expression\nX4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"
comparisons[["X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC"]][["contrast"]] <- c("condition","X4006.Osi_plusJH_10_day_or_DTC", "X4006.Osi_10_day_or_DTC")

### Comparison End ##


cts.file <- "/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/STAR_Gene_Counts.csv"
cts <- read.table(cts.file,sep=",",header=TRUE,stringsAsFactors=FALSE)

any(duplicated(rownames(cts$Gene_ID)))
rownames(cts) <- cts$Gene_ID
cts$Gene_ID <- NULL


universe.entrez <- na.omit(mapIds(org.Hs.eg.db, rownames(cts), 'ENTREZID', 'SYMBOL'))


comp.tables <- list()
for (comparison in names(comparisons)) {
  
    comp.name <- comparison

    cts.filtered <- cts[,comparisons[[comparison]][["columns"]]]
    
    coldata <- data.frame(condition=comparisons[[comparison]][["condition"]],
                          batch = comparisons[[comparison]][["batch"]])
    
    dds <- DESeqDataSetFromMatrix(countData = cts.filtered,
                                  colData = coldata,
                                  design= ~ condition)
    
    featureData <- data.frame(gene=rownames(cts.filtered))
    mcols(dds) <- DataFrame(mcols(dds),featureData)
    
    dds <- DESeq(dds)
    norm.counts <- counts(dds,normalized=TRUE)
    
    res <- results(dds,contrast=comparisons[[comparison]][["contrast"]])

    res2 <- cbind(norm.counts,res)
    
    pdf(file=paste(comparison,"_ma_plot.pdf",sep=""))
    plotMA(res,main=comparisons[[comparison]][["title"]])
    dev.off()
    
    resOrdered <- res[order(res$log2FoldChange),]
    resOrdered.df <- as.data.frame(resOrdered)
    resOrdered.df2 <- resOrdered.df[order(resOrdered.df$padj),]
    resOrdered.df2 <- resOrdered.df2[abs(resOrdered.df2$log2FoldChange) >= 1 & resOrdered.df2$padj <= 0.05 & (! is.na(resOrdered.df2$padj)),]
    padj.thresh <- min(0.05,resOrdered.df2[round(nrow(resOrdered.df2)/5),]$padj)

    up.genes <- rownames(resOrdered.df[which(resOrdered.df$padj <= 0.05 & resOrdered.df$log2FoldChange >= 1),])
    down.genes <- rownames(resOrdered.df[which(resOrdered.df$padj <= 0.05 & resOrdered.df$log2FoldChange <= -1),])
    
    gg.df <- data.frame(Log2FoldChange = resOrdered.df$log2FoldChange,
                        Log10PAdjValue = -log10(resOrdered.df$padj),
                        Gene=rownames(resOrdered.df),
                        DE="NO",
                        DE2="NO",
                        Labels="",
                        YAP1 = "", 
                        ERK = "")
    if (length(up.genes) > 0) {
      gg.df[gg.df$Gene %in% up.genes,]$DE <- "Up"
    }
    if (length(down.genes) > 0) {
      gg.df[gg.df$Gene %in% down.genes,]$DE <- "Down"
    }
    ggSig.df <- gg.df
    ggSig2.df <- gg.df
    gg.df[which(gg.df$Gene %in% c(up.genes, down.genes)),]$Labels <- gg.df[which(gg.df$Gene %in% c(up.genes, down.genes)),]$Gene
    pdf(file=paste(comparison,"_gg_vplot.pdf",sep=""),width=10)
    p <- ggplot(data=gg.df, aes(x=Log2FoldChange, y=Log10PAdjValue, col=DE, label=Labels)) +
      ylab("-Log10 Adj p-value")+
      geom_point() + theme_bw()+ theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                                       panel.background = element_blank())  +
      geom_text_repel(color="black",max.overlaps = 25) +
      geom_vline(xintercept=c(-1,1),col="black",linetype="dotted") +
      geom_hline(yintercept=-log10(0.05),col="black",linetype="dotted")
    if(length(up.genes) > 0 & length(down.genes) > 0) {
      p <- p + scale_color_manual(values=c(alpha("royalblue",0.5),alpha("gray",0.5),alpha("red3",0.5)))
    } else if (length(up.genes) > 0) {
      p <- p + scale_color_manual(values=c("gray","red3"))
    } else if (length(down.genes > 0 )) {
      p <- p + scale_color_manual(values=c("royalblue","gray"))
    } else {
      p <- p + scale_color_manual(values=c("gray"))
    }
    print(p)
    dev.off()
    
    # sig <- read.csv('/Users/nt788/Desktop/RNA_Seq/Bulk/20MSTO/signature_list120121.csv')
    # yap1.sig <- sig$YAP1_FINAL[which(sig$YAP1_FINAL != "")]
    # ggSig.df[ggSig.df$Gene %in% yap1.sig,]$Labels <- ggSig.df[ggSig.df$Gene %in% yap1.sig,]$Gene
    # ggSig.df$YAP1 <- ifelse(ggSig.df$Labels== "", "N/A", "YAP1")
    # pdf(file=paste(comparison,"_gg_YAP1_vplot.pdf",sep=""),width=10)
    # p <- ggplot(data=ggSig.df, aes(x=Log2FoldChange, y=Log10PAdjValue, col=YAP1, label=Labels)) +
    #   ylab("-Log10 Adj p-value")+
    #   geom_point() + theme_bw()+ theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
    #                                    panel.background = element_blank())  +
    #   geom_text_repel(color="black",max.overlaps = Inf) +
    #   geom_vline(xintercept=c(-1,1),col="black",linetype="dotted") +
    #   geom_hline(yintercept=-log10(0.05),col="black",linetype="dotted")
    #   p <- p + scale_color_manual(values = c("N/A" = alpha('lightgray', 0.1), "YAP1" = 'violet'))
    # 
    # 
    # print(p)
    # dev.off()
    # 

    
  
  # ggSig2.df[which(sig$ERK_FINAL %in% c(up.genes, down.genes)),]$Labels <- ggSig2.df[which(sig$ERK_FINAL  %in% c(up.genes, down.genes)),]$Gene
  # erk.sig <- sig$ERK_FINAL[which(sig$ERK_FINAL != "")]
  # ggSig2.df[ggSig2.df$Gene %in% erk.sig,]$Labels <- ggSig2.df[ggSig2.df$Gene %in% erk.sig,]$Gene
  # ggSig2.df$ERK <- ifelse(ggSig2.df$Labels== "", "N/A", "ERK")
  # pdf(file=paste(comparison,"_gg_ERK_vplot.pdf",sep=""),width=10)
  #   p <- ggplot(data=ggSig2.df, aes(x=Log2FoldChange, y=Log10PAdjValue, col=ERK, label=Labels)) +
  #     ylab("-Log10 Adj p-value")+
  #     geom_point() + theme_bw()+ theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
  #                                      panel.background = element_blank())  +
  #     geom_text_repel(color="black",max.overlaps = Inf) +
  #     geom_vline(xintercept=c(-1,1),col="black",linetype="dotted") +
  #     geom_hline(yintercept=-log10(0.05),col="black",linetype="dotted")
  #   p <- p + scale_color_manual( values = c("N/A" = alpha('lightgray', 0.1) , "ERK" = 'violet'))
  # 
  # 
  #   print(p)
  #   dev.off()
  #   
  

    write.table(resOrdered.df,file=paste(comparison,".txt",sep=""),sep="\t")
    write.table(up.genes,file=paste(comparison,"_up_genes.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
    write.table(down.genes,file=paste(comparison,"_down_genes.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
    comp.tables[[comparison]] <- as.data.frame(res2)

    if (length(up.genes) >= 10) {
        up.entrez <- na.omit(mapIds(org.Hs.eg.db, up.genes, 'ENTREZID', 'SYMBOL'))

        up.go.test <- enrichGO(up.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
        #up.go.test <- go.simplify(up.go.test)
        pdf(file=paste(comparison,"_up_go.pdf",sep=""),width=12,height=14)
        print(dotplot(up.go.test,showCategory=40))
        dev.off()
        write.csv(data.frame(up.go.test),(file = paste(comparison,"_up_go.csv",sep="")))

     }

    if (length(down.genes) >= 10) {
        down.entrez <- na.omit(mapIds(org.Hs.eg.db, down.genes, 'ENTREZID', 'SYMBOL'))

        down.go.test <- enrichGO(down.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
        #down.go.test <- go.simplify(down.go.test)
        pdf(file=paste(comparison,"_down_go.pdf",sep=""),width=12,height=14)
        print(dotplot(down.go.test,showCategory=40))
        dev.off()

        write.csv(data.frame(down.go.test),(file = paste(comparison,"_down_go.csv",sep="")))

    }

    log2.fds <- resOrdered.df$log2FoldChange
    names(log2.fds) <- rownames(resOrdered.df)
    gene.list <- na.omit(log2.fds)
    gene.list = sort(gene.list, decreasing = TRUE)

    gmt <- read.gmt("/Users/nt788/Desktop/RNA_Seq/Bulk/20MSTO/hallmark.gmt")
    pgmt <-read.gmt("/Users/nt788/Desktop/RNA_Seq/Bulk/20MSTO/c2.cp.v2022.1.Hs.symbols.gmt")
    egmt <- GSEA(gene.list, TERM2GENE=gmt,minGSSize =6,
                 pvalueCutoff = 0.05,pAdjustMethod = "fdr")
    pdf(paste0(comparison, "_gseaHallmark.pdf"),width=12,height=10)
    print(ggplot(egmt@result, aes(reorder(Description, NES), NES)) +
            geom_col(aes(fill=p.adjust)) +
            coord_flip() +
            geom_text(aes(label = signif(p.adjust,digits=2)), size=5, nudge_y  = -0.2, colour = "black") +
            labs(x="Pathway", y="Normalized Enrichment Score",
                 title=comparisons[[comparison]][["title"]]) +
            theme_minimal()+ theme(legend.text = element_text(colour="black", size=10),
                                   axis.text.y =element_text(colour="black", size=15),
                                   axis.title.y =element_text(colour="black", size=15),
                                   axis.text.x =element_text(colour="black", size=15),
                                   axis.title.x =element_text(colour="black", size=15))+
            scale_fill_gradient2(low = "red", high = "blue",mid = "purple",limits=c(0,1),midpoint = 0.5)
    )
    dev.off()

    # egmt <- GSEA(gene.list, TERM2GENE=pgmt,minGSSize =6,
    #              pvalueCutoff = 0.1,pAdjustMethod = "fdr")
    # pdf(paste0(comparison, "_gseaHumanPathways.pdf"),width=24, height=20)
    # print(top_n(egmt@result, n = -30, egmt@result$p.adjust) %>%
    #   ggplot(., aes(reorder(Description, NES), NES)) +
    #         geom_col(aes(fill=p.adjust)) +
    #         coord_flip() +
    #         geom_text(aes(label = signif(p.adjust,digits=2)), size=5, nudge_y  = -0.2, colour = "black") +
    #         labs(x="Pathway", y="Normalized Enrichment Score",
    #              title=comparisons[[comparison]][["title"]]) +
    #         theme_minimal()+ theme(legend.text = element_text(colour="black", size=10),
    #                                axis.text.y =element_text(colour="black", size=15),
    #                                axis.title.y =element_text(colour="black", size=15),
    #                                axis.text.x =element_text(colour="black", size=15),
    #                                axis.title.x =element_text(colour="black", size=15))+
    #         scale_fill_gradient2(low = "red", high = "blue",mid = "purple",limits=c(0,1),midpoint = 0.5)
    # )
    # dev.off()


    up.genes <- rownames(resOrdered.df[which(resOrdered.df$pvalue <= 0.05 & resOrdered.df$log2FoldChange >= 1),])
    down.genes <- rownames(resOrdered.df[which(resOrdered.df$pvalue <= 0.05 & resOrdered.df$log2FoldChange <= -1),])

    gg.df <- data.frame(Log2FoldChange = resOrdered.df$log2FoldChange,
                        Log10PValue = -log10(resOrdered.df$pvalue),
                        Gene=rownames(resOrdered.df),
                        DE="NO",
                        DE2="NO",
                        Labels="")
    if (length(up.genes) > 0) {
      gg.df[gg.df$Gene %in% up.genes,]$DE <- "Up"
    }
    if (length(down.genes) > 0) {
      gg.df[gg.df$Gene %in% down.genes,]$DE <- "Down"
    }
    gg.df[which(gg.df$Gene %in% c(up.genes, down.genes)),]$Labels <- gg.df[which(gg.df$Gene %in% c(up.genes, down.genes)),]$Gene

    pdf(file=paste(comparison,"_gg_pvalue_vplot.pdf",sep=""),width=10)
    p <- ggplot(data=gg.df, aes(x=Log2FoldChange, y=Log10PValue, col=DE, label=Labels)) +
      ylab("-Log10 p-value")+
      geom_point() + theme_bw()+ theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                                       panel.background = element_blank())  +
      geom_text_repel(color="black",max.overlaps = 25) +
      geom_vline(xintercept=c(-1,1),col="black",linetype="dotted") +
      geom_hline(yintercept=-log10(0.05),col="black",linetype="dotted")
    if(length(up.genes) > 0 & length(down.genes) > 0) {
      p <- p + scale_color_manual(values=c(alpha("royalblue",0.5),alpha("gray",0.5),alpha("red3",0.5)))
    } else if (length(up.genes) > 0) {
      p <- p + scale_color_manual(values=c("gray","red3"))
    } else if (length(down.genes > 0 )) {
      p <- p + scale_color_manual(values=c("royalblue","gray"))
    } else {
      p <- p + scale_color_manual(values=c("gray"))
    }
    print(p)
    dev.off()



     if (length(up.genes) >= 10) {
      up.entrez <- na.omit(mapIds(org.Hs.eg.db, up.genes, 'ENTREZID', 'SYMBOL'))

      up.go.test <- enrichGO(up.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
      #up.go.test <- go.simplify(up.go.test)
      pdf(file=paste(comparison,"_pvalue_up_go.pdf",sep=""),width=12,height=14)
      print(dotplot(up.go.test,showCategory=40))
      dev.off()
      ##write.csv(data.frame(up.go.test),(file = paste(comparison,"_up_go.csv",sep="")))

    }

    if (length(down.genes) >= 10) {
      down.entrez <- na.omit(mapIds(org.Hs.eg.db, down.genes, 'ENTREZID', 'SYMBOL'))

      down.go.test <- enrichGO(down.entrez, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05,ont="BP",universe=universe.entrez,maxGSSize=2000)
      #down.go.test <- go.simplify(down.go.test)
      pdf(file=paste(comparison,"_pvalue_down_go.pdf",sep=""),width=12,height=14)
      print(dotplot(down.go.test,showCategory=40))
      dev.off()

      ##write.csv(data.frame(down.go.test),(file = paste(comparison,"_down_go.csv",sep="")))

    }



}


