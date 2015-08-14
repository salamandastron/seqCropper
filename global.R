source("http://bioconductor.org/biocLite.R")
biocLite()
reqLibs <- c("ChIPpeakAnno","biomaRt","BSgenome","XML","mailR"
             #              ,"BSgenome.Hsapiens.NCBI.GRCh38",
             #              "BSgenome.Mmusculus.UCSC.mm9",
             #              "BSgenome.Mmusculus.UCSC.mm10",
             #              "BSgenome.Hsapiens.UCSC.hg18",
             #              "BSgenome.Hsapiens.UCSC.hg19"
)

cat("Checking status of required packages and reference genomes...\n")

if (length(intersect(reqLibs,installed.packages()[,1])) < 5) {
  cat("Installing genomes and required packages...\n")
  #cat("If this is the first time you are using me, it might take a while to load all the stuff. We are working with genomes after all...")
  aaa = installed.packages()
  biocLite(setdiff(reqLibs,as.character(aaa[,1])))
}

cat("Required packages ready.\n\nPipeline ready to run.\n")

library(ChIPpeakAnno)
library(biomaRt)
library(mailR)

# Provisional input options
# upstreamLen <- 0
# downstreamLen <- 200
# dbName <- 'ensembl'
# genome <- 'human-GRCh38'
# genomeVers <- 'GRCh38'
# # Currently supports:
# # transcription start site, exons, transcripts, cDNA, 5'UTR, 3'UTR
# featureType <- 'TSS'
# # Supports Ensembl gene IDs and HGNC symbols as identifiers
# inputFormat <- 'symbol'
# inputQuery <- "MYC"
exon.num <- 1
earliest.flag <- T
reactiveVal <- reactiveValues(currmart=useMart("ensembl", dataset = 'hsapiens_gene_ensembl'),
                              idtype="hgnc_symbol")
# genomeLookup <- c("BSgenome.Hsapiens.NCBI.GRCh38",
#                   "BSgenome.Mmusculus.UCSC.mm9",
#                   "BSgenome.Mmusculus.UCSC.mm10",
#                   "BSgenome.Hsapiens.UCSC.hg18",
#                   "BSgenome.Hsapiens.UCSC.hg19")
# names(genomeLookup) <- c("human-GRCh38","mouse-mm9","mouse-mm10","human-hg18","human-hg19")
# library(genomeLookup[genome])


# if (length(grep('human',genome,fixed=T)) > 0) {
#   ds <- 'hsapiens_gene_ensembl'
# } else if (length(grep('human',genome,fixed=T)) > 0) {
#   ds <- 'mmusculus_gene_ensembl'
# } else {
#   cat("Ooops I don't have data for that species yet...maybe check back later?\n")
# }
# 
# mart = useMart("ensembl", dataset = ds)

getExon <- function(gene,exon.num,mart,idtype) {
  gb <- getBM(attributes=c('ensembl_gene_id','ensembl_exon_id','rank',
                           "exon_chrom_start","exon_chrom_end","gene_exon","strand"),
              filters = idtype, values=gene, mart=mart,
              bmHeader=TRUE)
  #gb <- gb[grep("ENSE",gb$`Ensembl Exon ID`),]
#   strandSign <- as.character(gb$Strand[1])
#   if (strandSign == "1") {
#     seqPos <- gb$`Exon Chr Start (bp)`
#     exonOrder <- order(seqPos,decreasing = F)
#   } else {
#     seqPos <- gb$`Exon Chr End (bp)`
#     exonOrder <- order(seqPos,decreasing = T)
#   }
  exonSeq <- gb$`Exon sequences`[gb$`Exon Rank in Transcript`==exon.num]
}

getTSS <- function(gene,earliest=T,mart,idtype) {
  #cat(".")
  gb <- getBM(attributes=c('transcription_start_site','chromosome_name','strand'),
              filters=idtype, values=gene, mart=mart)
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

getPromoter <- function(gene,upflank,mart,idtype) {
  #cat("...")
  #   gb <- getBM(attributes=c('transcription_start_site','chromosome_name','strand'),
  #               filters="hgnc_symbol", values=gene, mart=mart)
  resSeq <- getSequence(id = gene,type=idtype,
                        mart=mart,seqType="gene_flank",upstream=upflank)
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
  
  #   gb <- getBM(attributes=c('ensembl_gene_id','ensembl_exon_id',
  #                            "exon_chrom_start","exon_chrom_end","gene_exon","strand"),
  #               filters = "hgnc_symbol", values=gene, mart=mart,
  #               bmHeader=TRUE)
}

cat("\n")