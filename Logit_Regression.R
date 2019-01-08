require(plyr)
require(dplyr)

setwd('/Users/Mel/Desktop/CALeDNA/deco_3')
metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
metadata$Zeta_4 <- as.factor(as.character(metadata$Zeta_4)) #convert Zeta_4 to factor
metadata$clust <- as.factor(as.character(metadata$clust)) #convert clust to factor
character_vars <- lapply(metadata, class) == "character" #convert all char columns to factor
metadata[, character_vars] <- lapply(metadata[, character_vars], as.factor)
int_vars <- lapply(metadata, class) == "integer" #convert all int columns to numeric
metadata[, int_vars] <- lapply(metadata[,int_vars], as.numeric)

#function to convert logit coeffs. to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#Get the groups for the type of environmental factors.
factorGroups <- read.table("Metadata_explanation.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Unique environmental factors are as follows:
uniqueFactors <- c("Location","Topology","Habitat","BioClim","Soil Properties","Vegetation") #, "Human Impact")

#Get all of the metagenomic data tables, clustered by family, within a primer set.
primerList <- c("16S","18S","CO1") #,"PITS","FITS",
totalLRM <- data.frame()
for(primer in primerList){
  primerFile <- paste("ppsubset",primer,".csv",sep="")
  primerSubset <- read.table(primerFile,header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, row.names = 1)
  primerSubset <- as.data.frame(sapply(primerSubset, as.factor))
  uniqueSpecies <- c(colnames(primerSubset))
  if (ncol(primerSubset >0)){
    for(uniqueFactor in uniqueFactors){
      factorSubset <- subset(factorGroups,Category==uniqueFactor)
      metadataSubset <- metadata
      rownames(metadataSubset) <- metadataSubset$MatchName
      metadataSubset <- metadataSubset[,colnames(metadataSubset) %in% factorSubset$`Column name`]
      #Force metadata factors into numeric or factor types.
      metadataSubset[sapply(metadataSubset, is.character)] <- lapply(metadataSubset[sapply(metadataSubset, is.character)], as.factor)
      metadataSubset[sapply(metadataSubset, is.integer)] <- lapply(metadataSubset[sapply(metadataSubset, is.integer)], as.numeric)
      indivFactors <- colnames(metadataSubset) #individual factors that make up uniqueFactor group i.e. bio1 through bio 19
      #Estimate log regression model
      for(uniqueSpecie in uniqueSpecies){
        mergeDF <- cbind(select(primerSubset, uniqueSpecie),metadataSubset)
        names(mergeDF)[names(mergeDF)==uniqueSpecie] <- "y"
        if(length(unique(mergeDF$y)) > 1){
          LRM <- glm(y ~ . - y, data = mergeDF, family = binomial("logit"),maxit=1000, na.action = na.omit)
          LRMfit <- as.data.frame(coef(summary(LRM)))
          #convert logit coefficient (Estimate) to probability of presence
          LRMfit$Probability_of_Presence <- logit2prob(LRMfit$Estimate)
          LRMfit$taxaID <- uniqueSpecie
          LRMfit$primer <- primer
          LRMfit$factor <- row.names(LRMfit)
          LRMfit <- subset(LRMfit,LRMfit$`Pr(>|z|)`<=0.05)
          if('(Intercept)' %in% row.names(LRMfit)){
            LRMfit <- LRMfit[-1,]
          }
        }
        totalLRM <- rbind(totalLRM, LRMfit)
      }
    }
  }
}

write.csv(totalLRM, file = "/Users/Mel/Desktop/CALeDNA/deco_3/totalLRM.csv", row.names = TRUE)
