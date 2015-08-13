library(shiny)
shinyUI(fluidPage(
  titlePanel("seqCropper"),
  
  sidebarLayout(
    sidebarPanel(
      #helpText("Batch retrieve your sequence of interest."),
      selectInput("var", 
                  label = "Choose a genome",
                  choices = c("Mouse - GRCm38/mm10", "Human - GRCh38/hg38"),
                  selected = "Human - GRCh38/hg38"),
      textInput("text", label = h3("Your favorite gene"), 
                value = "MYC"),
      helpText("For promoter sequences, enter desired length upstream of TSS:"),
      textInput("upWindow",label="Upstream sequence length",
                value = "1000"),
      actionButton("submit","Get promoter sequences"),
      helpText("For exon sequences, enter which exon you are looking for:"),
      textInput("exonNum",label="Exon number",
                value = "1"),
      actionButton("submitExon","Get exon sequences"),
      helpText("For gene position information:"),
      actionButton("submitGenePos","Get gene positions"),
      helpText("For UTR sequences:"),
      radioButtons("utr", "UTR selection:",
                   c("5' UTR" = "5",
                     "3' UTR" = "3")),
      actionButton("submitUTR","Get UTR sequences")
#       sliderInput("range", 
#                   label = "Range of interest:",
#                   min = 0, max = 100, value = c(0, 100))
    ),
    
    mainPanel(
      verbatimTextOutput("genome"),
      tabsetPanel(type = "tabs", 
                  tabPanel("Promoter", verbatimTextOutput("promoter")), 
                  tabPanel("Exon", verbatimTextOutput("exon")), 
                  tabPanel("Gene Position", dataTableOutput("genePosTable")),
                  tabPanel("UTR", verbatimTextOutput("utr"))
      )
      #verbatimTextOutput("text")
      #textOutput("text")
    )
  )
  
))