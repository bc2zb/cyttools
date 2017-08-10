panel <- read.delim("CD11bPosCyttoolsResults/FlowSOM/cyttoolsBlanksResults/panelFile.txt")

panel$EQBEADS[grep("EQBEADS", panel$desc)] <- 1
panel$Barcode[grep("Pd", panel$desc)] <- 1
panel$EventDiscrimination[grep("Ir", panel$desc)] <- 1
panel$Viability[grep("Pt", panel$desc)] <- 1

panel$descOriginal <- panel$desc

panel$desc <- gsub(".*\\_|-EQBEADS", "", panel$desc)

panel$metal <- gsub("Di", "", panel$name)
panel$mass <- gsub("[A-z]*", "", panel$metal)
panel$element <- gsub("[0-9]*", "", panel$metal)
panel$nameCheck <- paste0(panel$mass, panel$element)  

panel$Ignore[which(panel$nameCheck == panel$desc)] <- 1
panel$Ignore[which(panel$Barcode == 1)] <- 1
panel$Ignore[grep("EQBEADS|dist", panel$desc)] <- 1
panel$Ignore[is.na(panel$desc)] <- 1
panel$Ignore[which(panel$Viability == 1)] <- 1
panel$Ignore[which(panel$EventDiscrimination == 1)] <- 1

functionalMarkerIndex <- grep("p", panel$desc)

panel$Lineage[panel$Ignore == 0 & panel$Functional == 0] <- 1

metaData <- read.delim("CD11bPosCyttoolsResults/FlowSOM/cyttoolsBlanksResults/MetaDataFile.txt")

metaData$realFileName <- gsub("CD11b DBC_24HR 25NG HBO23.fcs", "", metaData$FileName)

metaData <- separate(metaData,
                     realFileName,
                     colnames(metaData)[c(5,4,2)],
                     sep = "\ ")

metaData <- metaData[,c(1,8,3,7,6)]
metaData$SampleID <- gsub("\\.fcs", "", metaData$SampleID)
metaData$TimePoint <- gsub(".*\\_", "", metaData$TimePoint)
metaData$SampleID[19:20] <- rep("BatchControl", 2)
metaData$Condition[19:20] <- NA
metaData$TimePoint[19:20] <- NA
