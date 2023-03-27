library(VennDiagram)
library(DESeq2)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(UpSetR)
library(plyr)   


## Step 1: Set genes to index -- DONE
regTable <- read.csv('/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/comparisonTable.csv', row.names='X0')
regTable$emtSig <- 0
regTable$integrin <-0

EMT <- read.csv('//Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/EMT/HallMark_EMT_Sig.csv')
EMT$emtSig <- 1
rownames(EMT) <-EMT$gene
for(i in rownames(EMT)){
  regTable[i, 'emtSig'] <- 1
}

integrin <- read.delim('/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/String_integrinMediatedSignalingPathway.tsv', sep = '\t')
integrin$integrin <-1
rownames(integrin) <- integrin$X.node
for(i in rownames(integrin)){
  regTable[i, 'integrin'] <- 1
}




## Step 2: Create 3 tables-- where Up is 1 other is 0, Down is 1 other is 0, & NA is 1 other is 0 -- DONE
comps <- colnames(regTable)

upReg <- regTable
upReg[upReg == "Up"] <- 1
upReg[upReg == "Down"] <- 0
upReg[is.na(upReg)] <- 0
upReg <- upReg %>% mutate_if(is.character, as.numeric)



downReg <- regTable
downReg[downReg == "Up"] <- 0
downReg[downReg == "Down"] <- 1
downReg[is.na(downReg)] <- 0
downReg <- downReg %>% mutate_if(is.character, as.numeric)

naReg <- regTable
naReg[naReg == "Up"] <- 0
naReg[naReg == "Down"] <- 0
naReg[is.na(naReg)] <- 1
naReg <- naReg %>% mutate_if(is.character, as.numeric)

#Step 3: Plot
pdf(file=paste("upRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(upReg, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("downRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(downReg, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("naRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(naReg, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()


## Step 4: Refine Plots

# 24 Hour
up24 <- upReg[,c("PC9.OsiplusJH_24_vs_PC9.Osi_24", "X4006.OsiplusJH_24_vs_X4006.Osi_24",'integrin')]
down24 <- downReg[,c("PC9.OsiplusJH_24_vs_PC9.Osi_24", "X4006.OsiplusJH_24_vs_X4006.Osi_24", 'integrin')]
na24 <- naReg[,c("PC9.OsiplusJH_24_vs_PC9.Osi_24", "X4006.OsiplusJH_24_vs_X4006.Osi_24")]
comps <- c("PC9.OsiplusJH_24_vs_PC9.Osi_24", "X4006.OsiplusJH_24_vs_X4006.Osi_24",'integrin')

pdf(file=paste("24hrUpRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(up24, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("24hrDownRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(down24, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("24hrNaRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(na24, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()


# 72 Hour
up72 <- upReg[,c("PC9.OsiplusJH_72_vs_PC9.Osi_72", "X4006.OsiplusJH_72_vs_X4006.Osi_72",  'integrin')]
down72 <- downReg[,c("PC9.OsiplusJH_72_vs_PC9.Osi_72", "X4006.OsiplusJH_72_vs_X4006.Osi_72",  'integrin')]
na72 <- naReg[,c("PC9.OsiplusJH_72_vs_PC9.Osi_72", "X4006.OsiplusJH_72_vs_X4006.Osi_72")]
comps <- c("PC9.OsiplusJH_72_vs_PC9.Osi_72", "X4006.OsiplusJH_72_vs_X4006.Osi_72", 'integrin')

pdf(file=paste("72hrUpRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(up72, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("72hrDownRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(down72, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("72hrNaRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(na72, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

# 10 Day
up10 <- upReg[,c("PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC", "X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC",  'integrin')]
down10 <- downReg[,c("PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC", "X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC",  'integrin')]
na10 <- naReg[,c("PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC", "X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC")]
comps <- c("PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC", "X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC", 'integrin')

pdf(file=paste("10dayUpRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(up10, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("10dayDownRegulatedUpset_Integrin.pdf",sep=""),width=15,height=12)
print(upset(down10, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

pdf(file=paste("10dayNaRegulatedUpset.pdf",sep=""),width=15,height=12)
print(upset(na10, sets = comps, empty.intersections = "on", order.by = "freq"))
dev.off()

## Step 5: Pull gene lists
list24 <- regTable[,c("PC9.OsiplusJH_24_vs_PC9.Osi_24", "X4006.OsiplusJH_24_vs_X4006.Osi_24", 'integrin')]
up24List <- subset(list24, PC9.OsiplusJH_24_vs_PC9.Osi_24 == 'Up' & X4006.OsiplusJH_24_vs_X4006.Osi_24 == 'Up' & integrin == 1)
down24List <- subset(list24, PC9.OsiplusJH_24_vs_PC9.Osi_24 == 'Down' & X4006.OsiplusJH_24_vs_X4006.Osi_24 == 'Down'&  integrin == 1)
write.csv(up24List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/24HrUpRegulatedOverlap_Integrin.csv')
write.csv(down24List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/24HrDownRegulatedOverlap_Integrin.csv')


list72 <- regTable[,c("PC9.OsiplusJH_72_vs_PC9.Osi_72", "X4006.OsiplusJH_72_vs_X4006.Osi_72", 'integrin')]
up72List <- subset(list72, PC9.OsiplusJH_72_vs_PC9.Osi_72 == 'Up' & X4006.OsiplusJH_72_vs_X4006.Osi_72 == 'Up' &  integrin == 1)
down72List <- subset(list72, PC9.OsiplusJH_72_vs_PC9.Osi_72 == 'Down' & X4006.OsiplusJH_72_vs_X4006.Osi_72 == 'Down'&  integrin == 1)
write.csv(up72List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/72HrUpRegulatedOverlap_Integrin.csv')
write.csv(down72List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/72HrDownRegulatedOverlap_Integrin.csv')


list10 <- regTable[,c("PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC", "X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC",  'integrin')]
up10List <- subset(list10, PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC == 'Up' & X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC == 'Up'& integrin == 1)
down10List <- subset(list10, PC9.Osi_plusJH_10_day_or_DTC_vs_PC9.Osi_10_day_or_DTC == 'Down' & X4006.Osi_plusJH_10_day_or_DTC_vs_X4006.Osi_10_day_or_DTC == 'Down'& integrin == 1)
write.csv(up10List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/10DayUpRegulatedOverlap_Integrin.csv')
write.csv(down10List, '/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/10DayDownRegulatedOverlap_Integrin.csv')

# emtIntegrin <- regTable[,c('integrin', 'emtSig')]
# emtIntegrinList <- subset(emtIntegrin, integrin == 1 & emtSig == 1)
# write.csv(emtIntegrinList,'/Users/nt788/Desktop/RNA_Seq/Bulk/OsiPlusJH/overlappingGenes/EMT_IntegrinOverlap.csv' )

