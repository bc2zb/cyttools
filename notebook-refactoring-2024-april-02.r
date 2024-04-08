#' ---
#' title: "cyttools refactoring"
#' author: "Brian Capaldo"
#' date: "4/2/2024"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: hide
#' ---

#+ message = FALSE, warning = FALSE

#### environment creation ####

#' # Prototyping environment
#'
#' Will be creating a conda env using mamba to get things up and running as
#' quickly as possible. Reviewing code, definitely some inefficiencies to be
#' addressed, but right now, focus is on getting the project running again.
#' Extensive refactoring will be taking place. `cyttoolsFunctions.R` contains
#' all the libraries of packages required.

#+ message = FALSE, warning = FALSE

package_list <- c("radian",
                  "r-base",
                  "r-essentials",
                  "r-devtools",
                  "r-biocmanager",
                  "bioconductor-limma",
                  "bioconductor-ncdfflow",
                  "bioconductor-opencyto",
                  "bioconductor-biocparallel",
                  "bioconductor-edger",
                  "bioconductor-flowsom",
                  "bioconductor-flowcore")

create_string <- paste(c("mamba create --name cyttools_env",
    package_list),
    collapse = " ")

#### test for libraries ####

libraryList <- c("flowCore",
                 "limma",
                 "edgeR",
                 "FlowSOM",
                 "ncdfFlow",
                 "openCyto",
                 "tidyverse")

lapply(libraryList, require, quietly = T, character.only = TRUE)

#### add in docopt ####

install.packages("docopt")

#### export yml ####

export_string <- "conda env export > cyttools.yml"

#### test flowVS code ####

source("cyttoolsFunctions.R")


dir <- "/data/capaldobj/cidc-cyttools-test-bed/"
file <- list.files(dir,
                   pattern = '.fcs$',
                   full = TRUE)

#### use panel blank code as template ####
library(flowVS)
flowSet <- read.ncdfFlowSet(file)

panel <- as.data.frame(pData(parameters(flowSet[[1]])))

cofactors_est <- flowVS::estParamFlowVS(flowSet, colnames(flowSet)[-c(1:2)])
dev.off()

targets <- targets %>% mutate(Lineage = if_else(desc %in% c("Time",
                                                 "Event_length",
                                                 "Viability",
                                                 "DNA") | desc == name,
                                        0,
                                        1))
