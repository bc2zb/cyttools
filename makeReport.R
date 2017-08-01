#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
makeReport.R (-h | --help | --version)
makeReport.R DIR

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory where cyttools.args.Rdata file is located
" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# diff abundance analysis

diffAbndncFile <- paste0(args$DIFFDIR, "nodeDifferentialAbundanceTable.txt")
diffAbndncResults <- read.delim(diffAbndncFile)
diffAbndncResults <- diffAbndncResults[complete.cases(diffAbndncResults),]

# diff expr analysis

diffExprFile <- paste0(args$DIFFDIR, "nodeDifferentialExpressionTable.txt")
diffExprResults <- read.delim(diffExprFile)
diffExprResults <- diffExprResults[complete.cases(diffExprResults),]

# grab results that are sig in one, other or both

popsIdsFile <- paste0(args$CLUSTERDIR, "PhenoCodes.txt")
popsIds <- read.delim(popsIdsFile)

rm(diffAbndncFile)
rm(diffExprFile)
rm(popsIdsFile)

diffAbndncResults$Mapping <- recoderFunc(diffAbndncResults$Mapping, popsIds$Names, popsIds$PhenoCodes)

sigAbndncNodes <- diffAbndncResults$Mapping[diffAbndncResults$adj.P.Val < 0.01]
sigDiffNodes <- diffExprResults[diffExprResults$adj.P.Val < 0.01 ,c(8:9)]

sigNodes <- c(sigAbndncNodes, sigDiffNodes$Mapping) %>% unique()
sigNodes <- factor(sigNodes, levels = popsIds$PhenoCodes)
sigNodes <- sigNodes[order(sigNodes, decreasing = T)]

numNodesToPlot <- 512
NodesToPlot <- c(1:numNodesToPlot)

sigAbndncResults <- diffAbndncResults[diffAbndncResults$Mapping %in% sigNodes[NodesToPlot],]
sigExprResults <- diffExprResults[diffExprResults$Mapping %in% sigNodes[NodesToPlot],]

sigAbndncResults$Phenotypes <- factor(recoderFunc(sigAbndncResults$Mapping,
                                                  popsIds$PhenoCodes,
                                                  popsIds$Names),
                                      levels = popsIds$Names)

ggplot(sigAbndncResults, aes(Source, Phenotypes, fill = logFC)) +
  geom_tile() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "ReportWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
