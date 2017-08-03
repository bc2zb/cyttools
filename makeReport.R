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

load("~/BTSync/FLOW_CORE/Rivera-Nieves/ColonPBMCPaper/cyttoolsReport/cyttools.args.Rdata")# 
RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

# ##########################################################################
# ############################ R code goes here ############################
# ##########################################################################

# diff abundance analysis

diffAbndncFile <- paste0(args$DIFFDIR, "nodeDifferentialAbundanceTable.txt")
diffAbndncResults <- read.delim(diffAbndncFile)
diffAbndncResults <- diffAbndncResults[complete.cases(diffAbndncResults),]

# diff expr analysis

diffExprFile <- paste0(args$DIFFDIR, "nodeDifferentialExpressionTable.txt")
diffExprResults <- read.delim(diffExprFile)
diffExprResults <- diffExprResults[complete.cases(diffExprResults),]

# grab results that are sig in one, other or both DOESN'T EXIST FOR FLOWSOM RESULTS

#popsIdsFile <- paste0(args$CLUSTERDIR, "PhenoCodes.txt")
#popsIds <- read.delim(popsIdsFile)

# panel <- read.delim(args$PANEL)
# targets <- read.delim(args$METADATA)
# 
# rm(diffAbndncFile)
# rm(diffExprFile)
# rm(popsIdsFile)
# 
# diffAbndncResults$Mapping <- recoderFunc(diffAbndncResults$Mapping, popsIds$Names, popsIds$PhenoCodes)
# diffAbndncResults$Metal <- rep("Abundance", nrow(diffAbndncResults))
# 
# diffResults <- rbind(diffExprResults, diffAbndncResults)
# 
# diffResults$Mapping <- factor(diffResults$Mapping, levels = popsIds$PhenoCodes)
# diffResults$Metal <- recoderFunc(diffResults$Metal,
#                                  panel$name,
#                                  gsub(".*_|-EQBEADS", "", panel$desc))
# 
# diffResults <- separate(diffResults,
#                         Source,
#                         c("Hold", "Vary"),
#                         "\\.",
#                         remove = F)
# 
# # p1 <- ggplot(diffResults,
# #        aes(Metal,
# #            Mapping,
# #            fill = -log10(adj.P.Val)
# #            )
# #        ) +
# #   geom_tile() + 
# #   facet_grid(Hold ~ Vary) +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1),
# #         axis.text.y = element_blank()
# #         ) + 
# #   scale_fill_gradient2(low = "grey0",
# #                        high = "red",
# #                        mid = "white",
# #                        midpoint = -log10(0.10))
# # 
# # tiff(file = paste0(args$OUT, "allResultsHeatMap.tiff"),
# #     width = 10.24,
# #     height = 7.68,
# #     res = 600,
# #     units = "in")
# # print(p1)
# # dev.off()
# # 
# # pdf(file = paste0(args$OUT, "allComparisonsHeatMaps.pdf"),
# #      width = 10.24,
# #      height = 7.68)
# # 
# # for(i in 1:length(unique(diffResults$Source))){
# #   sigSourceResults <- diffResults[diffResults$Source == unique(diffResults$Source)[i],]
# #   sigSourceResults <- sigSourceResults[sigSourceResults$adj.P.Val < 0.10,]
# #   p2 <- ggplot(sigSourceResults[sigSourceResults$Metal != "CD11b",],
# #                aes(Metal,
# #                    Mapping,
# #                    fill = logFC
# #                )
# #   ) +
# #     geom_tile() + 
# #     theme(axis.text.x = element_text(angle = 90, hjust = 1),
# #           axis.text.y = element_blank()
# #     ) + 
# #     scale_fill_gradient2(low = "blue",
# #                          high = "red",
# #                          mid = "white",
# #                          na.value = "grey",
# #                          midpoint = 0) + 
# #     ggtitle(as.character(unique(diffResults$Source)[i]))
# #   print(p2)
# # }
# # dev.off()
# 
# popsIds$filledPhenotypes <- sapply(popsIds$PhenoCodes, function(x){
#   x <- as.character(x)
#   len <- nchar(x)
#   if(len < 9){
#     x <- paste0(paste(rep(0, (9 - len)),
#                       collapse = ""),
#                 x,
#                 collapse = "")    
#   }
#   return(x)
# })
# 
# popsIds <- extract(popsIds,
#                    filledPhenotypes,
#                    popsIds$Names[2:19] %>% gsub(".*\\_|\\+|\\-|-EQBEADS", "", .) %>% unique(),
#                    regex = '(.)(.)(.)(.)(.)(.)(.)(.)(.)',
#                    remove = F)
# 
# ParentPops <- popsIds[grep("CD11b|CD14-|CD16", popsIds$Names[1:19]),]
# 
# popsIds$NumMarkers <- apply(popsIds[,c(4:ncol(popsIds))], 1, function(x){x %>%  gsub("2", "1", .) %>% as.numeric() %>% sum()})
# popsIds$PhenoCodes <- factor(popsIds$PhenoCodes, levels = popsIds$PhenoCodes)
# 
# diffResults <- diffResults  %>%
#                left_join(popsIds, by = c("Mapping" = "PhenoCodes"))
# 
# parentMarkers <- c("CD11b", "CD16", "CD14")
# childMarkers <- c("CD184", "CD38", "CD142", "CD86", "CD43", "CD197")
# 
# subPopDiffResults <- dplyr::filter(diffResults,
#                              NumMarkers == 4 &
#                              Metal == "Abundance") %>%
#                      gather(ChildMarker,
#                             ExprScore, 
#                             CD184, CD38, CD142, CD86, CD43, CD197) %>%
#                      dplyr::filter(ExprScore != 0 &
#                                    CD11b != 0 &
#                                    CD16 != 0 &
#                                    CD14 != 0)
# 
# subPopDiffResults$Names <- gsub("-EQBEADS|\\d*\\w*\\_", "", subPopDiffResults$Names)
# subPopDiffResults$ChildNames <- gsub("CD11b\\-|CD14\\-|CD16\\-|CD11b\\+|CD14\\+|CD16\\+", "", subPopDiffResults$Names)
# subPopDiffResults$ParentNames <- apply(subPopDiffResults, 1, function(x){
#   x <- gsub(x[names(x) %in% "ChildNames"], "", x[names(x) %in% "Names"], fixed = T)
#   return(x)
# })
# 
# subPopDiffResults$sigAnnotation <- subPopDiffResults$adj.P.Val
# subPopDiffResults$sigAnnotation[subPopDiffResults$adj.P.Val < 0.01] <- "***"
# subPopDiffResults$sigAnnotation[subPopDiffResults$adj.P.Val < 0.05 &
#                                 subPopDiffResults$adj.P.Val >= 0.01] <- "**"
# subPopDiffResults$sigAnnotation[subPopDiffResults$adj.P.Val < 0.1 &
#                                   subPopDiffResults$adj.P.Val >= 0.05] <- "**"
# subPopDiffResults$sigAnnotation[subPopDiffResults$adj.P.Val >= 0.1] <- ""
# plot <- ggplot(subPopDiffResults[subPopDiffResults$Source %in% levels(subPopDiffResults$Source)[c(1:3,5:6,8:9)],],
#        aes(ChildNames,
#            logFC,
#            # ,
#            # fill = logFC,
#            fill = sigAnnotation
#        )
# ) +
#   geom_bar(stat = "identity") + 
#   #geom_text() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4.5),
#         legend.position = "none",
#         strip.text.y = element_text(angle = 0),
#         strip.text.x = element_text(size = 4.5)) +
#   xlab("") + 
#   ylab("") +
#   facet_grid(Source ~ ParentNames) +
#   scale_fill_manual(values = c("grey0", "black", "red"))
# 
# plotFile <- paste(RESULTS_DIR, "SummaryAbundancePlot.pdf", sep = "")
# pdf(file = plotFile,
#     width = 10.24,
#     height = 7.68)
# print(plot)
# dev.off()

##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "ReportWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
