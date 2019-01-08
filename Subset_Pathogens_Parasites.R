require(dplyr)

setwd('/Users/Mel/Desktop/CALeDNA/deco_3')

pList = read.csv("pList.csv") #pathogen and parasite list

#read in OTU tables
species16S = read.csv("asv_deco_dedup_16S.csv")
species18S = read.csv("asv_deco_dedup_18S.csv")
speciesCO1 = read.csv("asv_deco_dedup_CO1.csv")
speciesFITS = read.csv("asv_deco_dedup_FITS.csv")
speciesPITS = read.csv("asv_deco_dedup_PITS.csv")

#initialize empty data frames
subset16S <- data.frame()
subset18S <- data.frame()
subsetCO1 <- data.frame()
subsetFITS <- data.frame()
subsetPITS <- data.frame()

#filter OTU by partial match 
for(row in 1:nrow(pList)){
  match <- species16S %>% dplyr::filter(grepl(as.character(pList[row,1]), sum.taxonomy)) #16S
  # print(paste(as.character(pList[row,1]),nrow(match)))
  if(nrow(match) >= 1){
    subset16S <- rbind(match, subset16S)}
  
  match <- species18S %>% dplyr::filter(grepl(as.character(pList[row,1]), sum.taxonomy)) #18S
  # print(paste(as.character(pList[row,1]),nrow(match)))
  if(nrow(match) >= 1){
    subset18S <- rbind(match, subset18S)}
  
  match <- speciesCO1 %>% dplyr::filter(grepl(as.character(pList[row,1]), sum.taxonomy)) #CO1
  # print(paste(as.character(pList[row,1]),nrow(match)))
  if(nrow(match) >= 1){
    subsetCO1 <- rbind(match, subsetCO1)}

  match <- speciesFITS %>% dplyr::filter(grepl(as.character(pList[row,1]), sum.taxonomy)) #FITS
  # print(paste(as.character(pList[row,1]),nrow(match)))
  if(nrow(match) >= 1){
    subsetFITS <- rbind(match, subsetFITS)}
  
  match <- speciesPITS %>% dplyr::filter(grepl(as.character(pList[row,1]), sum.taxonomy)) #PITS
  # print(paste(as.character(pList[row,1]),nrow(match)))
  if(nrow(match) >= 1){
    subsetPITS <- rbind(match, subsetPITS)}
    }

# convert subset to presence/absence counts i.e. all numbers >= 1 become 1
subset16S[subset16S >= 1] <- 1
subset18S[subset18S >= 1] <- 1
subsetCO1[subsetCO1 >= 1] <- 1
subsetFITS[subsetFITS >= 1] <- 1
subsetPITS[subsetPITS >= 1] <- 1

#transform subset OTUs
subset16S = setNames(data.frame(t(subset16S[,-1])), subset16S[,1])
subset18S = setNames(data.frame(t(subset18S[,-1])), subset18S[,1])
subsetCO1 = setNames(data.frame(t(subsetCO1[,-1])), subsetCO1[,1])
subsetFITS = setNames(data.frame(t(subsetFITS[,-1])), subsetFITS[,1])
subsetPITS = setNames(data.frame(t(subsetPITS[,-1])), subsetPITS[,1])

#rename headers in subset OTUs to species only
names(subset16S) <- sub(".* ", "", names(subset16S))
names(subset18S) <- sub(".* ", "", names(subset18S))
names(subsetCO1) <- sub(".* ", "", names(subsetCO1))
names(subsetFITS) <- sub(".* ", "", names(subsetFITS))
names(subsetPITS) <- sub(".* ", "", names(subsetPITS))

#save subset OTU primer tables as csv
write.csv(subset16S, file = "/Users/Mel/Desktop/CALeDNA/deco_3/ppsubset16S.csv", row.names = TRUE)
write.csv(subset18S, file = "/Users/Mel/Desktop/CALeDNA/deco_3/ppsubset18S.csv", row.names = TRUE)
write.csv(subsetCO1, file = "/Users/Mel/Desktop/CALeDNA/deco_3/ppsubsetCO1.csv", row.names = TRUE)
write.csv(subsetFITS, file = "/Users/Mel/Desktop/CALeDNA/deco_3/ppsubsetFITS.csv", row.names = TRUE)
write.csv(subsetPITS, file = "/Users/Mel/Desktop/CALeDNA/deco_3/ppsubsetPITS.csv", row.names = TRUE)
