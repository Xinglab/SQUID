args <- commandArgs(trailingOnly = TRUE)
### using intron and exon data as input
library(DEXSeq)
##agrs[1] = intron_data
##args[2] = sampleData
##args[3] = alternativeCountData

intron = read.table(args[1],row.names = 1)
rownames(intron) = paste(intron[,1], rownames(intron),sep=":")
intron = intron[,-1]
countData = intron
featureID = gsub(".*:","",rownames(countData))
groupID = gsub(":.*","",rownames(countData))
groupID = gsub(",.*","",groupID)

sampleData = read.table(args[2])
sampleData = as.data.frame(sampleData)
alternativeCountData = read.table(args[3])
colnames(alternativeCountData) = colnames(countData)
alternativeCountData  = as.matrix(alternativeCountData )
dxd = DEXSeqDataSet(countData, sampleData, design= ~ sample + exon + condition:exon , featureID, groupID, featureRanges=NULL, transcripts=NULL, alternativeCountData)

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )

dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
write.table(dxr1,args[length(args)],quote = FALSE, sep="\t")
