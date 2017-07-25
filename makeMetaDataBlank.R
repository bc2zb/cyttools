#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
makeMetaDataBlank.R (-h | --help | --version)
makeMetaDataBlank.R DIR

Description:   This script is makes meta data table blanks for use with cyttools.
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory where fcs files are located
" -> doc


args <- docopt(doc)

RESULTS_DIR <- args$DIR

cat("\nLoading arguments from", args$DIR, "\n")

load(paste(RESULTS_DIR, "cyttools.args.Rdata", sep = ""))

library(tidyverse)

# USER SHOULD BE PASSING THIS IN AS A SEPARATE FILE
dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

md <- data.frame(FileName = gsub(paste(args$DIR, "/", sep = "|"), "", file))
md <- data.frame(FileName = md$FileName,
                 SampleID = vector(length = length(md$FileName)),
                 Group = vector(length = length(md$FileName)),
                 Condition = vector(length = length(md$FileName)),
                 TimePoint = vector(length = length(md$FileName)))
      
mdFile <- paste(RESULTS_DIR, "MetaDataFile.txt", sep = "")

write.table(md, file = mdFile, sep = "\t", quote = F, row.names = F)
