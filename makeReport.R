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

panel <- read.delim(args$PANEL)
targets <- read.delim(args$METADATA)

rm(diffAbndncFile)
rm(diffExprFile)
rm(popsIdsFile)

diffAbndncResults$Mapping <- recoderFunc(diffAbndncResults$Mapping, popsIds$Names, popsIds$PhenoCodes)
diffAbndncResults$Metal <- rep("Abundance", nrow(diffAbndncResults))

diffResults <- rbind(diffExprResults, diffAbndncResults)

diffResults$Mapping <- factor(diffResults$Mapping, levels = popsIds$PhenoCodes)
diffResults$Metal <- recoderFunc(diffResults$Metal,
                                 panel$name,
                                 gsub(".*_|-EQBEADS", "", panel$desc))

diffResults <- separate(diffResults,
                        Source,
                        c("Hold", "Vary"),
                        "\\.",
                        remove = F)

p1 <- ggplot(diffResults,
       aes(Metal,
           Mapping,
           fill = -log10(adj.P.Val)
           )
       ) +
  geom_tile() + 
  facet_grid(Hold ~ Vary) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank()
        ) + 
  scale_fill_gradient2(low = "grey0",
                       high = "red",
                       mid = "white",
                       midpoint = -log10(0.10))

tiff(file = paste0(args$OUT, "allResultsHeatMap.tiff"),
    width = 10.24,
    height = 7.68,
    res = 600,
    units = "in")
print(p1)
dev.off()

pdf(file = paste0(args$OUT, "allComparisonsHeatMaps.pdf"),
     width = 10.24,
     height = 7.68)

for(i in 1:length(unique(diffResults$Source))){
  sigSourceResults <- diffResults[diffResults$Source == unique(diffResults$Source)[i],]
  sigSourceResults <- sigSourceResults[sigSourceResults$adj.P.Val < 0.10,]
  p2 <- ggplot(sigSourceResults[sigSourceResults$Metal != "CD11b",],
               aes(Metal,
                   Mapping,
                   fill = logFC
               )
  ) +
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_blank()
    ) + 
    scale_fill_gradient2(low = "blue",
                         high = "red",
                         mid = "white",
                         na.value = "grey",
                         midpoint = 0) + 
    ggtitle(as.character(unique(diffResults$Source)[i]))
  print(p2)
}
dev.off()

popsIds$filledPhenotypes <- sapply(popsIds$PhenoCodes, function(x){
  len <- nchar(x)
  if(len < 9){
    x <- paste0(paste(rep(0, (9 - len)),
                      collapse = ""),
                x,
                collapse = "")    
  }
  return(x)
})

codeMap <- extract(popsIds,
                   filledPhenotypes,
                   popsIds$Names[2:19] %>% gsub(".*\\_|\\+|\\-|-EQBEADS", "", .) %>% unique(),
                   regex = '(.)(.)(.)(.)(.)(.)(.)(.)(.)',
                   remove = F)

ParentPops <- codeMap[grep("CD11b|CD14-|CD16", popsIds$Names[1:19]),]

codeVals <- ParentPops[,c(4:ncol(codeMap))] %>% as.matrix() 
class(codeVals) <- "numeric"

ParentPopSeach <- colSums(codeVals) %>% paste0(., collapse = "") %>% gsub("0", ".", .) %>% gsub("3", "[1-2]", .)

diffResults$FullMap <- recoderFunc(diffResults$Mapping, popsIds$PhenoCodes, popsIds$filledPhenotypes)

subPopDiffResults <- diffResults[grep(ParentPopSeach, diffResults$FullMap),]
subPopDiffResults$ParentPopPhenoCode <- subPopDiffResults$FullMap
  
str_sub(subPopDiffResults$ParentPopPhenoCode[1], start = 2, end = 6) <- 0

# sigAbndncNodes <- diffAbndncResults$Mapping[diffAbndncResults$adj.P.Val < 0.01]
# sigDiffNodes <- diffExprResults[diffExprResults$adj.P.Val < 0.01 ,c(8:9)]
# 
# sigNodes <- c(sigAbndncNodes, sigDiffNodes$Mapping) %>% unique()
# sigNodes <- factor(sigNodes, levels = popsIds$PhenoCodes)
# sigNodes <- sigNodes[order(sigNodes, decreasing = T)]
# 
# numNodesToPlot <- 512
# NodesToPlot <- c(1:numNodesToPlot)
# 
# sigAbndncResults <- diffAbndncResults[diffAbndncResults$Mapping %in% sigNodes[NodesToPlot],]
# sigExprResults <- diffExprResults[diffExprResults$Mapping %in% sigNodes[NodesToPlot],]
# 
# sigAbndncResults$Phenotypes <- factor(recoderFunc(sigAbndncResults$Mapping,
#                                                   popsIds$PhenoCodes,
#                                                   popsIds$Names),
#                                       levels = popsIds$Names)

##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "ReportWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
