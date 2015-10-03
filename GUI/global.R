source("http://bioconductor.org/biocLite.R")
#biocLite()
reqLibs <- c("ChIPpeakAnno","biomaRt","BSgenome","XML","ggplot2",
             "TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene","ggbio","biovizBase"
             #              ,"BSgenome.Hsapiens.NCBI.GRCh38",
             #              "BSgenome.Mmusculus.UCSC.mm9",
             #              "BSgenome.Mmusculus.UCSC.mm10",
             #              "BSgenome.Hsapiens.UCSC.hg18",
             #              "BSgenome.Hsapiens.UCSC.hg19"
)

cat("Checking status of required packages and reference genomes...\n")

if (length(intersect(reqLibs,installed.packages()[,1])) < length(reqLibs)) {
  cat("Installing genomes and required packages...\n")
  #cat("If this is the first time you are using me, it might take a while to load all the stuff. We are working with genomes after all...")
  aaa = installed.packages()
  biocLite(setdiff(reqLibs,as.character(aaa[,1])))
}

if (!("mailR" %in% installed.packages()[,1])) {
  install.packages("mailR")
}

cat("Required packages ready.\n\nPipeline ready to run.\n")

library(ChIPpeakAnno)
library(biomaRt)
library(mailR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(ggbio)
library(biovizBase)

exon.num <- 1
earliest.flag <- T
reactiveVal <- reactiveValues(currmart=useMart("ensembl", dataset = 'hsapiens_gene_ensembl'),
                              idtype="hgnc_symbol",species = "Human")

drawGene <- function(gene,mart,sp,idtype) {
  genePos <- makeGRangesFromDataFrame(cbind.data.frame(getGenePos(gene,mart,idtype)))
  if (sp == "Human") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  gr.txdb <- crunch(txdb, which = genePos)
  colnames(values(gr.txdb))[4] <- "model"
  grl <- split(gr.txdb, gr.txdb$tx_id)
  names(grl) <- sapply(c(1:length(table(gr.txdb$tx_id))),function(x){
    paste0(gene,"_transcript_",as.character(x))
  })
  autoplot(grl, aes(type = model))
}

getExon <- function(gene,exon.num,mart,idtype) {
  gb <- getBM(attributes=c('ensembl_gene_id','ensembl_exon_id','rank',
                           "exon_chrom_start","exon_chrom_end","gene_exon","strand"),
              filters = idtype, values=gene, mart=mart,
              bmHeader=TRUE)
  exonSeq <- gb$`Exon sequences`[gb$`Exon Rank in Transcript`==exon.num]
  exonSeq
}

getGenePos <- function(gene,mart,idtype) {
  gb <- getBM(attributes=c('start_position','end_position',
                           'chromosome_name','strand'),
              filters=idtype, values=gene, mart=mart)
  strandSign <- as.character(gb$strand[1])
  if (strandSign == "1") {
    chr <- paste0("chr",gb$chromosome_name[1])
    startPos <- gb$start_position
    endPos <- gb$end_position
    strnd <- "+"
  } else {
    chr <- paste0("chr",gb$chromosome_name[1])
    startPos <- gb$end_position
    endPos <- gb$start_position
    strnd <- "-"
  }
  list(chr=chr,start=startPos,end=endPos,strand=strnd)
}

getTSS <- function(gene,earliest=T,mart,idtype) {
  #cat(".")
  gb <- getBM(attributes=c('transcription_start_site',
                           'chromosome_name','strand'),
              filters=idtype, values=gene, mart=mart)
  strandSign <- as.character(gb$strand[1])
  TSSList <- gb$transcription_start_site
  #TESList <- gb$transcription_end_site
  TSSOrder <- order(TSSList,decreasing=F)
  if (strandSign == "1") {
    if (earliest) {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList[TSSOrder[1]]
      #endPos <- TESList[TSSOrder[1]]
    } else {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList
      #endPos <- TESList
    }
  } else {
    if (earliest) {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList[length(TSSOrder)]
      #endPos <- TESList[length(TSSOrder)]
    } else {
      chr <- gb$chromosome_name[1]
      startPos <- TSSList
      #endPos <- TESList
    }
  }
  list(TSS=startPos,chr=chr,strand=strandSign)
}

getPromoter <- function(gene,upflank,mart,idtype) {
  #cat("...")
  #   gb <- getBM(attributes=c('transcription_start_site','chromosome_name','strand'),
  #               filters="hgnc_symbol", values=gene, mart=mart)
  resSeq <- getSequence(id = gene,type=idtype,
                        mart=mart,seqType="gene_flank",upstream=upflank)
  resSeq$gene_flank
}


getUTR <- function(gene,ori,mart,idtype) {
  if (ori == '5') {
    tag <- '5utr'
    resSeq <- getSequence(id = gene,type=idtype,mart=mart,seqType="5utr")
  } else {
    tag <- '3utr'
    resSeq <- getSequence(id = gene,type=idtype,mart=mart,seqType="3utr")
  }
  toKick <- grep("Sequence unavailable",resSeq[,1])
  resSeq <- resSeq[-toKick,]
  if (length(resSeq[[tag]])) {
    resSeq[[tag]]
  } else {
    "No specified UTR found for this gene..."
  }
}

cat("\n")