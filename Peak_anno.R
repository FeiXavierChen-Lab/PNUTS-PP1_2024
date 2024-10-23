library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

txdb_ucsc <- makeTxDbFromGFF(
  file = "/share/home/Blueberry/reference/annotation/ucsc/hg19.ncbiRefSeq.gtf",
  format = "gtf",
  organism = "Homo sapiens",
  taxonomyId = NA,
  dataSource = "ucscgenomes",
  metadata = NULL
)

files <- list.files("/chenfeilab/Coconut_srv1/PNUTS/Analysis/ChIPseq/ChIP_INTAC_SSB1/ssb1_used/", pattern = "Peak", full.names = T)
files_1 <- c(
  "/chenfeilab/Coconut_srv1/PNUTS/Analysis/ChIPseq/01_peakCompare_INTACandPNUTS/IP-FLAG_dTAG-3h_PNUTS-fkbp-DLD1_201208_compare_peaks.final.narrowPeak",
  "/chenfeilab/Coconut_srv1/PNUTS/Analysis/ChIPseq/01_peakCompare_INTACandPNUTS/IP-PNUTS_dTAG-3h_PNUTS-fkbp-DLD1_201208_compare_peaks.final.narrowPeak",
  files
)

peakAnnoList <- lapply(files_1, annotatePeak, TxDb = txdb_ucsc,
                       tssRegion = c(-1000, 1000), verbose = FALSE)
names(peakAnnoList) <- c("FLAG", "PNUTS", 
                         "INTS3_1", "INTS3_2", "INTS5_1", 
                         "INTS5_2", "SSB_1", "SSB_2")

plotAnnoBar(peakAnnoList)
peakAnnoList_1 <- peakAnnoList

for (i in names(peakAnnoList_1)) {
  peakAnnoList_1[[i]]@annoStat <- data.frame(
    Feature = c("Promoter", "5' UTR", "3' UTR", "Intragenic", "Intergenic"),
    Frequency = c(
      peakAnnoList_1[[i]]@annoStat[1, 2],
      peakAnnoList_1[[i]]@annoStat[2, 2],
      peakAnnoList_1[[i]]@annoStat[3, 2],
      sum(peakAnnoList_1[[i]]@annoStat[4:7, 2]),
      sum(peakAnnoList_1[[i]]@annoStat[8:9, 2])
    )
  )
  
  peakAnnoList_1[[i]]@annoStat$Feature <- factor(
    peakAnnoList_1[[i]]@annoStat$Feature,
    levels = c("Promoter", "5' UTR", "3' UTR", "Intragenic", "Intergenic")
  )
}

pdf("/chenfeilab/Coconut_srv1/PNUTS/Analysis/ChIPseq/01_peakCompare_INTACandPNUTS/Peak_anno_barchart.pdf")
plotAnnoBar(peakAnnoList_1) + theme_classic()
dev.off()
