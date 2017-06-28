library(flowType)

# Eventually ask for directory
# dir <- "~/BTSync/FLOW_CORE/ERICKSON_CYTOF_B_CELL_ALLERGY_CYTOF_DATA/"
# file <- list.files(dir ,pattern='.fcs$', full=TRUE)
# fileToSampleID <- read.delim("~/BTSync/FLOW_CORE/ERICKSON_DATA/FileNameToSampleIDTable.txt")
# 
# LABELS <- data.frame(x = gsub("^.*\\/\\/", "", file),
#                      fileName = file)
# 
# LABELS <- base::merge(LABELS, fileToSampleID, by.x = 1, by.y = 1)
# LABELS$Group <- colsplit(LABELS$sample.ID, "\ ", c("Group", "Uselss"))[,1]
# 
# gsub("^.*\\/\\/", "", file) %in% fileToSampleID$x
# flowSet <- read.flowSet(file)
# 
# targets <- cbind(colnames(flowSet), unlist(parameters(flowSet[[1]])$name), unlist(parameters(flowSet[[1]])$desc))
# 
# flowTypeTargets <- data.frame(PropMarkers = targets[,3],
#                               MFIMarkers = targets[,3],
#                               MarkerNames = targets[,3]
#                               )
# 
# 
# flowTypeTargets$MarkerNames <- gsub("\\ \\(v\\)", "", flowTypeTargets$MarkerNames)
# flowTypeTargets$MarkerNames <- gsub("-[A-Z][a-z][0-9]*|-EQ", "", flowTypeTargets$MarkerNames)
# 
# markersOfInterest <- c("IgM", "IgD", "HLA-DR", "CD27", "CXCR5", "CD24")
# variantMarkers <- c("CD196", "HLA-DR", "IgD", "CD45", "CD25", "CD43", "CD279-FITC", "CD38", "CD24", "CXCR5", "IgM", "CXCR4")
# 
# flowSet.trans <- transFlowVS(flowSet, 
#                              targets[flowTypeTargets$PropMarkers,1],
#                              rep(5, length(targets[flowTypeTargets$PropMarkers,1])))
# 
# load(file = "~/BTSync/FLOW_CORE/ERICKSON/DATA/flowTypeResultsAllEventsAllFilesMeatAllergyProjectUpdate.Rdata")
# 
# 
# # MFIs=Res@MFIs;
# # Proportions=Res@CellFreqs;
# # Proportions <- Proportions / max(Proportions)
# # names(Proportions) <- unlist(lapply(Res@PhenoCodes,
# #                                     function(x){return(decodePhenotype(x,
# #                                                                        as.character(Res@MarkerNames),
# #                                                                        Res@PartitionsPerMarker))}))
# # rownames(MFIs)=names(Proportions)
# # index=order(Proportions,decreasing=TRUE)[1:20]
# # bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
# # text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0,
# #      labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
# # axis(2);
# # axis(1, at=bp, labels=FALSE);
# # title(xlab='Phenotype Names', ylab='Cell Proportion')
# 
# phenotype.names=unlist(lapply(ResList[[1]]@PhenoCodes,function(x){
#   return(decodePhenotype(x,
#                          as.character(flowTypeTargets$MarkerNames[flowTypeTargets$MarkerNames %in% variantMarkers == T]),
#                          ResList[[1]]@PartitionsPerMarker))}))
# names(ResList[[1]]@PhenoCodes)=phenotype.names
# 
# all.proportions <- matrix(0,length(ResList[[1]]@CellFreqs),length(file))
# for (i in 1:length(ResList))
#   all.proportions[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]
# 
# Pvals <- vector(length = dim(all.proportions)[1]);
# EffectSize <- vector(length = dim(all.proportions)[1]);
# for (i in 1:dim(all.proportions)[1]){
#   
#      ##If all of the cell proportions are 1 (i.e., the phenotype
#      ##with no gates) the p-value is 1.
#     if (length(which(all.proportions[i,]!=1))==0){
#        Pvals[i]=1;
#        EffectSize[i]=0;
#        next;
#        }
#    temp=perm.t.test(all.proportions[i, LABELS$Group=="AD"],
#                 all.proportions[i, LABELS$Group=="HD"],
#                 statistic = "mean")
#  Pvals[i] <- temp$p.value
#  EffectSize[i] <- abs(temp$statistic)
#  }
# Pvals[is.nan(Pvals)]=1
# names(Pvals)=phenotype.names
# 
# save.image(file = "~/BTSync/FLOW_CORE/ERICKSON/DATA/FlowTypeStatsTestingUpdate.Rdata")
# 
# # selected <- which(p.adjust(Pvals) < 0.01)
# # print(names(selected))
# 
# # res<-RchyOptimyx(ResList[[1]]@PhenoCodes, -log10(Pvals),
# #                  ResList[[1]]@PhenoCodes[selected], factorial(6),FALSE)