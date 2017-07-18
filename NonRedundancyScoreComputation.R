library(flowCore)
library(cytofCore)
library(flowVS)
library(flowType)
library(tidyverse)
library(reshape2)
source("~/BTSync/MethylationGeneExpressionProfilePaper/SCRIPTS/expressionSetFunctions.R")


dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory
flowSet <- read.flowSet(file, transformation = F) # reads in files as flowSet, required for flowType

targets <- read.delim(args$PANEL)

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

flowSet.trans <- transFlowVS(flowSet, 
                             as.character(targets$name[which(targets$Lineage == 1 |
                                                               targets$Functional == 1)]),
                             rep(5, length(targets$name[which(targets$Lineage == 1 |
                                                                targets$Functional == 1)])))


# force user to pass in file WHITEBOARD!
# targets <- pData(parameters(flowSet[[1]]))
# write.table(targets, file = "~/BTSync/FLOW_CORE/LOUGHRAN/DATA/panel.txt", quote = F, row.names = F, sep = "\t")



## Define a function that calculates the NRS per sample
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(flowSet.trans[,colnames(flowSet.trans) %in% lineage_markers], NRS, use.exprs = TRUE)

nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

# make metadata, should be submitted with job

md <- unique(ggdf$sample_id) %>% 
  colsplit(.,
           pattern = "\ ",
           names = c("Sample", "Condition"))
md$FileName <- unique(ggdf$sample_id)
md$Condition <- gsub("\\.fcs", "", md$Condition)
md$Condition[c(22:25)] <- rep(md$Condition[5], 4)

ggdf <- melt(nrs_sample, id.var = "sample_id",
             value.name = "nrs", variable.name = "antigen")

ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
ggdf$condition <- recoderFunc(ggdf$sample_id, md$FileName, md$Condition)
ggdf$condition <- factor(ggdf$condition, levels = c("Norm", "WT", "Y640F", "D661Y", "NK"))

ggplot(ggdf, aes(x = antigen, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
