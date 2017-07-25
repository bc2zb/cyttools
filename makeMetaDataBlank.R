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
