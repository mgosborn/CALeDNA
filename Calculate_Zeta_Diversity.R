require(zetadiv)
require(plyr)

setwd('/Users/Mel/Desktop/to_levi/')

#read in 16S OTU table 
species16S = read.csv("asv_deco_dedup_16S.csv")

#read in CALeDNA site metadata
metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Calculate and add column for number of samples per location. This will be used as a cut off for calculating zeta diversity
locationFrequency <- as.data.frame(table(metadata$loc))
colnames(locationFrequency) <- c("loc","Freq")
metadata <- plyr::join(metadata,locationFrequency,by=c("loc"))

#zetaLocation <- data.frame()
zetaNum=4

#Get OTU names
taxa <- species16S[,1]
#Transpose the OTU table and merge it with the sample metadata
Tspecies16S <- as.data.frame(t(species16S[,-c(1)]))
#convert to presence/absence counts i.e. all numbers >= 1 become 1
Tspecies16Smerge[Tspecies16Smerge >= 1] <- 1
colnames(Tspecies16S) <- taxa
Tspecies16S$MatchName <- rownames(Tspecies16S)
Tspecies16Smerge <-plyr::join(Tspecies16S,metadata,by=c("MatchName"))
#Filter out locations with less than 4 samples.
Tspecies16Smerge <- subset(Tspecies16Smerge,Freq>=zetaNum)
#Tspecies16Smerge <- subset(Tspecies16Smerge,locationFrequency>=zetaNum)





