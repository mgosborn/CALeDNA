require(dplyr)
require(tidyr)

setwd('/Users/Mel/Desktop/CALeDNA/deco_3')

#read in parasite/pathogen list, OTU tables & family p/a table
pList = read.csv("pList.csv") #parasite/pathogen list
species16S = read.csv("asv_deco_dedup_16S.csv")
species18S = read.csv("asv_deco_dedup_18S.csv")
speciesCO1 = read.csv("asv_deco_dedup_CO1.csv")
FamilyPA = read.table("CALeDNAFamilyPA.txt", header = TRUE, sep = "\t", as.is = T, skip = 0,fill=TRUE,check.names=FALSE)
#combine data tables
megaFile <- rbind(species16S,species18S,speciesCO1)

#initialize empty data frames
megaFileSub <-data.frame()
subsetFPA <- data.frame()

#subset by partial match to pList
primerList <-c("16S","18S","CO1")
for(primer in primerList){
  for(row in 1:nrow(pList)){
    match <- megaFile %>% dplyr::filter(grepl(as.character(pList[row,1]),sum.taxonomy))
    if(nrow(match) >=1){
      megaFileSub <- rbind(match,megaFileSub)
    }
  }
}

fList <- data.frame(do.call('rbind',strsplit(as.character(megaFileSub$sum.taxonomy),';',fixed=TRUE)))
fList <- data.frame(fList[,4])
fList <- data.frame(fList[!duplicated(fList),])
names(fList) <- "Family" 
fList <- data.frame(fList[!grepl("NA",fList$Family),])

#filter Family PA table by partial match and save as subsetFPA
for(row in 1:nrow(fList)){
  match <- FamilyPA %>% dplyr::filter(grepl(as.character(fList[row,1]), Family))
  if(nrow(match) >= 1){
    subsetFPA <- rbind(match, subsetFPA)}
}

#transform table
subsetFPA = setNames(data.frame(t(subsetFPA[,-1])), subsetFPA[,1])

#rename headers in subset OTUs to family only
names(subsetFPA) <- sub(".*;", "", names(subsetFPA))

#save subset table as csv
write.csv(subsetFPA, file = "/Users/Mel/Desktop/CALeDNA/deco_3/FPAsubset.csv", row.names = TRUE)
