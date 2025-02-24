source("R/data_processing.R")
source("R/file_handlers.R")
source("R/plotting.R")

source("R/modules/tcga_module.R")
source("R/modules/pcawg_module.R")
source("R/modules/driver_mutation_module.R")

# Libraries
library(reticulate)
library(ggplot2)
library(dplyr)
library(purrr)
library(plotly)
library(readr)
library(shinyjs)
library(DT)
library(tidyr)

# Global variables
global_TCGA_PCAWG <- "dataSource"

# Load data
sample_PCAWG <- read.csv(file.path(global_TCGA_PCAWG, "CliPP_PCAWG.tsv"), sep='\t')
sample_TCGA <- read.csv(file.path(global_TCGA_PCAWG, "CliPP_TCGA.tsv"), sep='\t')
driver_mutation_data <- read.delim("./driver_mut_data.tsv")

# Data preprocessing
sample_PCAWG$dataset <- "PCAWG"
sample_TCGA$dataset <- "TCGA"
sample_data <- rbind(sample_PCAWG, sample_TCGA)

# Fix cancer type names
driver_mutation_data$cancer <- gsub("RCCC", "KIRC", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("COREAD", "CRC", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("AML", "LAML", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("CM", "SKCM", driver_mutation_data$cancer)