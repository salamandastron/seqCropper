####### Please un-comment the following lines if you do not have ######
####### the specified packages installed. #############################
# source("http://bioconductor.org/biocLite.R")
# biocLite("RCytoscape")
# install.packages("igraph")
# install.packages("shiny")
# install.packages("devtools")
# devtools::install_github("shiny-incubator", "rstudio")
library(shiny)
#library(shinyIncubator)
setwd("~/Tools/genomics/seqCropper/")
runApp("GUI")