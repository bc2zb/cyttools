#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
makeMetaDataBlank.R (-h | --help | --version)
makeMetaDataBlank.R DIR

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory
" -> doc


args <- docopt(doc)

RESULTS_DIR <- args$DIR

cat("\nLoading arguments from", args$DIR, "\n")

load(paste(RESULTS_DIR, "cyttools.args.Rdata", sep = ""))


# USER SHOULD BE PASSING THIS IN AS A SEPARATE FILE
md <- data.frame(FileName = colnames(countTable))
md$FileName <- gsub("/Volumes/som-fccf/FC_McNamaraCA/McSkimming/170602 Capaldo HK exp stim cells fcs files/Myeloid export//", "", md$FileName)
md <- separate(md,
               FileName,
               into = c("Subset", "TimePoint", "Condition", "SampleID"),
               sep = "\ ",
               remove = F
)
md$SampleID <- gsub("\\.fcs", "", md$SampleID)

workspaceFile <- paste(RESULTS_DIR, "makeMetaDataBlankWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
