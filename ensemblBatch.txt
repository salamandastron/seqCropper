rm(list=ls())
source("http://bioconductor.org/biocLite.R")
biocLite()

library(shiny)

# Prototype pipeline for pulling out sequences of interest
# using the biomaRt R interface and the Shiny GUI.

reqLibs <- c("ChIPpeakAnno","biomaRt","BSgenome"
#              ,"BSgenome.Hsapiens.NCBI.GRCh38",
#              "BSgenome.Mmusculus.UCSC.mm9",
#              "BSgenome.Mmusculus.UCSC.mm10",
#              "BSgenome.Hsapiens.UCSC.hg18",
#              "BSgenome.Hsapiens.UCSC.hg19"
             )

cat("Checking status of required packages and reference genomes...\n")

if (length(intersect(reqLibs,installed.packages()[,1])) < 3) {
  cat("Installing genomes and required packages...\n")
  #cat("If this is the first time you are using me, it might take a while to load all the stuff. We are working with genomes after all...")
  biocLite(setdiff(reqLibs,as.character(installed.packages()$Package)))
}

cat("Required packages ready.\n\nPipeline ready to run.\n")

library(ChIPpeakAnno)
library(biomaRt)

# Provisional input options
upstreamLen <- 0
downstreamLen <- 200
dbName <- 'ensembl'
genome <- 'human-GRCh38'
genomeVers <- 'GRCh38'
# Currently supports:
# transcription start site, exons, transcripts, cDNA, 5'UTR, 3'UTR
featureType <- 'TSS'
# Supports Ensembl gene IDs and HGNC symbols as identifiers
inputFormat <- 'symbol'
inputQuery <- "MYC"
exon.num <- 1
earliest.flag <- T

# genomeLookup <- c("BSgenome.Hsapiens.NCBI.GRCh38",
#                   "BSgenome.Mmusculus.UCSC.mm9",
#                   "BSgenome.Mmusculus.UCSC.mm10",
#                   "BSgenome.Hsapiens.UCSC.hg18",
#                   "BSgenome.Hsapiens.UCSC.hg19")
# names(genomeLookup) <- c("human-GRCh38","mouse-mm9","mouse-mm10","human-hg18","human-hg19")
# library(genomeLookup[genome])


if (length(grep('human',genome,fixed=T)) > 0) {
  ds <- 'hsapiens_gene_ensembl'
} else if (length(grep('human',genome,fixed=T)) > 0) {
  ds <- 'mmusculus_gene_ensembl'
} else {
  cat("Ooops I don't have data for that species yet...maybe check back later?\n")
}

getExon <- function(gene,exon.num) {
  gb <- getBM(attributes=c('ensembl_gene_id','ensembl_exon_id',
                           "exon_chrom_start","exon_chrom_end","gene_exon","strand"),
              filters = "hgnc_symbol", values=gene, mart=mart,
              bmHeader=TRUE)
  strandSign <- as.character(gb$Strand[1])
  if (strandSign == "1") {
    seqPos <- gb$`Exon Chr Start (bp)`
    exonOrder <- order(seqPos,decreasing = F)
  } else {
    seqPos <- gb$`Exon Chr End (bp)`
    exonOrder <- order(seqPos,decreasing = T)
  }
  exonSeq <- gb$`Exon sequences`[exonOrder[exon.num]]
}

getTSS <- function(gene,earliest=T) {
  gb <- getBM(attributes=c('transcription_start_site','chromosome_name','strand'),
        filters="hgnc_symbol", values=gene, mart=mart)
#   gb <- getBM(attributes=c('gene_flank','start_position','end_position','chromosome_name',
#                      'strand','ensembl_gene_id','transcription_start_site'),
#               filters=c('hgnc_symbol','upstream_flank'),
#         values=list(gene,Upstream=1), mart=mart, checkFilters=FALSE)
  strandSign <- as.character(gb$strand[1])
  TSSList <- gb$transcription_start_site
  #TSSList <- as.numeric(strsplit(TSSList,";",fixed=T)[[1]])
  TSSOrder <- order(TSSList,decreasing=F)
  if (strandSign == "1") {
    if (earliest) {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList[TSSOrder[1]]
    } else {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList
    }
  } else {
    if (earliest) {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList[length(TSSOrder)]
    } else {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList
    }
  }
  list(TSS=startPos,chr=chr,strand=strandSign)
}

getPromoter <- function(gene,upflank) {
#   gb <- getBM(attributes=c('transcription_start_site','chromosome_name','strand'),
#               filters="hgnc_symbol", values=gene, mart=mart)
  resSeq <- getSequence(id = gene,type="hgnc_symbol",mart=mart,seqType="gene_flank",upstream=upflank)
  resSeq$gene_flank
#   gb <- getBM(attributes=c('gene_flank','start_position','end_position','chromosome_name',
#                        'strand','transcription_start_site'),
#                 filters=c('hgnc_symbol','upstream_flank'),
#           values=list(gene,Upstream=upflank), mart=mart, checkFilters=FALSE)
#   strandSign <- as.character(gb$strand[1])
#   TSSList <- gb$transcription_start_site
#   #TSSList <- as.numeric(strsplit(TSSList,";",fixed=T)[[1]])
#   TSSOrder <- order(TSSList,decreasing=F)
#   if (strandSign == "1") {
#     if (earliest) {
#       chr <- gb$chromosome_name[1]
#       startPos <- TSSList[TSSOrder[1]]
#     } else {
#       chr <- gb$chromosome_name[1]
#       startPos <- TSSList
#     }
#   } else {
#     if (earliest) {
#       chr <- gb$chromosome_name[1]
#       startPos <- TSSList[length(TSSOrder)]
#     } else {
#       chr <- gb$chromosome_name[1]
#       startPos <- TSSList
#     }
#   }
#   list(TSS=startPos,chr=chr,strand=strandSign)
}


getUTR <- function(gene,ori) {
  if (ori == '5') {
    resSeq <- getSequence(id = gene,type="hgnc_symbol",mart=mart,seqType="5utr")
  } else {
    resSeq <- getSequence(id = gene,type="hgnc_symbol",mart=mart,seqType="3utr")
  }
  toKick <- grep("Sequence unavailable",resSeq[,1])
  resSeq <- resSeq[-toKick,]
  resSeq
#   gb <- getBM(attributes=c('ensembl_gene_id','ensembl_exon_id',
#                            "exon_chrom_start","exon_chrom_end","gene_exon","strand"),
#               filters = "hgnc_symbol", values=gene, mart=mart,
#               bmHeader=TRUE)
}

cat("\n")

if (featureType == 'exon') {
  mart = useMart("ensembl", dataset = ds)
  resList <- list()
  for (g in gene) {
    resList[[g]] <- getExon(g,exon.num)
    cat(".")
  }
  cat("\n")
} else if (featureType == 'TSS') {
  mart = useMart("ensembl", dataset = ds)
  resList <- list()
  for (g in gene) {
    resList[[g]] <- getTSS(g,earliest=earliest.flag)
    cat(".")
  }
  cat("\n")
} else if (featureType == 'promoter') {
  mart = useMart("ensembl", dataset = ds)
  resList <- list()
  for (g in gene) {
    resList[[g]] <- getPromoter(g,upstreamLen)
    cat(".")
  }
  cat("\n")
} else if (featureType == "5'UTR") {
  mart = useMart("ensembl", dataset = ds)
  resList <- list()
  for (g in gene) {
    resList[[g]] <- getTSS(g,"5")
    cat(".")
  }
  cat("\n")
} else if (featureType == "3'UTR") {
  mart = useMart("ensembl", dataset = ds)
  resList <- list()
  for (g in gene) {
    resList[[g]] <- getTSS(g,"3")
    cat(".")
  }
  cat("\n")
} else {
  cat("Feature currently not supported...maybe check back later?\n")
}



getLDS(attributes = c("hgnc_symbol","chromosome_name", "start_position"),
       filters = "hgnc_symbol", values = "TP53",mart = mart)

annot.tss <- getAnnotation(mart,featureType="TSS")
annot.df <- as.data.frame(annot.tss)
genes <- read.delim("~/Downloads/results (5).txt",header=T)
geneIDs <- as.character(genes$Ensembl.Gene.ID)
annot.df.ROI <- annot.df[match(geneIDs,annot.df$names),]
#annot.df.ROI <- annot.df.ROI[-toDel,]
regionList <- annot.df.ROI[,c(1,2,3,5,6)]
# regionList$space <- sapply(as.character(regionList$space),function(x){
#   paste0("chr",x)
# })
toDel <- grep("CHR_",annot.df.ROI$space,fixed=T)
toDel <- c(toDel,grep("LRG_",annot.df.ROI$space,fixed=T))
regionList <- regionList[-toDel,]
write.table(regionList,"~/Desktop/rl.txt",row.names=F,sep="\t",quote=F)
regionList <- read.table("~/Desktop/rl.txt",header=T)
windowSize <- 200
seqList <- apply((regionList),1,function(x){
  chr <- as.character(x[1])
  #centerPos <- as.numeric(x[2])
  strandSign <- as.numeric(x[4]) # Do not know the strand actually...
  if (strandSign == 1) {
    strnd <- "+"
    startPos <- as.numeric(x[2])
    endPos <- startPos + windowSize - 1
  } else {
    strnd <- "-"
    endPos <- as.numeric(x[3])
    startPos <- endPos - windowSize + 1
  }
  #c(chr,as.character(startPos),as.character(endPos),strnd)
  cat(".")
  getSeq(BSgenome.Hsapiens.NCBI.GRCh38,chr,
         startPos,endPos,
         strand=strnd,as.character=T)
})
rgn.df <- t(seqList)
write.table(rgn.df,"~/Desktop/rgn.txt",row.names=F,sep="\t",quote=F)





