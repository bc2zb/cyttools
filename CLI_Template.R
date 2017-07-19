#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
CLI_Template.R (-h | --help | --version)
CLI_Template.R DIR

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

workspaceFile <- paste(RESULTS_DIR, "TemplateWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
