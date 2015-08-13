shinyServer(
  function(input, output) {
    
    setGenome <- reactive({
      if (input$var=="Human - GRCh38/hg38") {
        reactiveVal$currmart <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
        reactiveVal$idtype <- "hgnc_symbol"
        #getPromoter(input$text,as.numeric(input$upWindow),mart,idtype)
      } else if (input$var=="Mouse - GRCm38/mm10") {
        reactiveVal$currmart <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
        reactiveVal$idtype <- "mgi_symbol"
        #getPromoter(input$text,as.numeric(input$upWindow),mart,idtype)
      } else {
        message("Species not supported...")
      }
      input$var
    })
    geneInput.prom <- eventReactive({
      input$submit
      }, 
      {
        getPromoter(input$text,as.numeric(input$upWindow),reactiveVal$currmart,reactiveVal$idtype)
    })
    geneInput.exon <- eventReactive({
      input$submitExon
      }, 
      {
        getExon(input$text,as.numeric(input$exonNum),reactiveVal$currmart,reactiveVal$idtype)
      })
    geneInput.UTR <- eventReactive({
      input$submitUTR
      }, 
      {
        getUTR(input$text,(input$utr),reactiveVal$currmart,reactiveVal$idtype)
      })
    geneInput.TSS <- eventReactive({
      input$submitGenePos
    }, 
    {
      getTSS(input$text,mart=reactiveVal$currmart,idtype=reactiveVal$idtype)
    })
#     output$text1 <- renderText({ 
#       paste("You have selected", input$var)
#     })
#     output$ds <- renderText({    
#       setGenome()
#     })
    output$genome <- renderText({
      paste(setGenome(),collapse="\n")
    })
    output$promoter <- renderText({    
      paste(geneInput.prom(),collapse="\n")
    })
    output$exon <- renderText({    
      paste(geneInput.exon(),collapse="\n")
    })
    output$utr <- renderText({    
      paste(geneInput.UTR(),collapse="\n")
    })
    output$genePosTable <- renderDataTable(({
      cbind.data.frame(geneInput.TSS())
    }))
  }
)
