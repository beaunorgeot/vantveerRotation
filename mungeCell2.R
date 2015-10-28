setwd("/Users/beaunorgeot/vantveerRotation")
#Headers are in the 3rd row
mainTable <- read.csv("/Users/beaunorgeot/vantveerRotation/SuppTab1.csv", skip = 2)

library(dplyr)
#remove info I know I won't use
mainTab <- mainTable %>% select(-c(ExomeSeq.availability,RNASeq.availability,Exon.array.availability,U133A.availability,
                                   SNP6.availability,Methylation.availability, RPPA.availability))

#only take observations actually used in drug analysis
mainTab <- mainTab %>% filter(Included.in.drug.analyses == 1)
mainTab <- mainTab %>% select(-c(Included.in.drug.analyses))

#quick mean check to determine if values from meanGI50 are the means from this table
mainTab %>% select(X17.AAG) %>% summarise(mean(X17.AAG, na.rm=T)) # 7.096
# value from meanGI50 is 7.0354:: these values are not the same
# how many NAs were there?
summary(mainTab$X17.AAG) #19, thats 27% (19/70)
#try another drug
mainTab %>% select(Sigma.AKT1.2.inhibitor) %>% summarise(mean(Sigma.AKT1.2.inhibitor,na.rm=T)) #5.437
# value from meanGI50 is 5.46. Weird

#Did they replace NAs w/0??
zero4na <- mainTab %>% mutate(X17.AAG = ifelse(is.na(X17.AAG) ==T,0,X17.AAG)) #mean = 5.17, no they didn't replace na w/0
#what about replacing w/the mean?
#zero4na <- mainTab %>% mutate(X17.AAG = ifelse(is.na(X17.AAG) ==T,mean(X17.AAG,na.rm=T),X17.AAG))
# obviously, the above returns the same mean as w/o this transformation
# they also didn't replace NAs w/the median
zero4na <- mainTab %>% mutate(X17.AAG = ifelse(is.na(X17.AAG) ==T,median(X17.AAG),X17.AAG)) #mean = 7.096

#Add variables for sensitivity/res
drugs <- mainTab %>% select(-c(Cell.line,Transcriptional.subtype,Transcriptional.subtype...ERBB2.status))
cellDesciptions <- mainTab %>% select(Cell.line,Transcriptional.subtype,Transcriptional.subtype...ERBB2.status)
features <- names(drugs)

# median impute
for (f in features) {
  drugs[,f] <- ifelse(is.na(drugs[,f]), median(drugs[,f], na.rm= TRUE), drugs[,f])
}

drugs1 <- drugs
# This works!
#drugs1 <- drugs1 %>% mutate(Sens_X17.AAG = ifelse(X17.AAG < mean(X17.AAG,na.rm=T),1,0))

# Create Sensitivity DummyVar using median
for(f in features){
  dummy.var <- paste("Sens", f, sep='_')
  drugs1[, dummy.var] <- NA
  drugs1[,dummy.var] <- ifelse(drugs1[,f] <= median(drugs1[,f],na.rm=T),1,0) 
}

# Create Resistant DummyVar using mean
for(f in features){
  dummy.var <- paste("Res", f, sep='_')
  drugs1[, dummy.var] <- NA
  drugs1[,dummy.var] <- ifelse(drugs1[,f] >= median(drugs1[,f],na.rm=T),1,0) 
}

bidrugs <- drugs1[,91:270]
# just create a count summary
sumBidrugs <- bidrugs %>% summarise_each(funs(sum))
#sumBidrugs <- as.data.frame(t(sumBidrugs))
#sumBidrugs <- as.data.frame(t(sumBidrugs), row.names = NULL)
#sumBidrugs2 <- cbind(ResponseOfDrug = rownames(sumBidrugs), numObservations = sumBidrugs$V1)

#create 2 x J table tracking senstive and resistant counts for each drug
drugTable <- data.frame(matrix(NA, nrow = 2, ncol = 90))
colnames(drugTable) <- names(drugs1[,1:90])
rownames(drugTable) <- c("Sensitive", "Resistant")

for (i in 1:90){
  drugTable[1,i] <- sumBidrugs[1,i]
  drugTable[2,i] <- sumBidrugs[1,i+90]
}


save(drugTable, file = "drugContigencyTable.Rdata")


# EVERYTHING BELOW HERE BELONGS IN EXPLORE
avgCounts <- rowMeans(drugTable)
#Sensitive Resistant 
#46.32222  47.58889
library(dplyr)
#how many drugs have either an especially high/low ratio of sens/res?

responseProps <- as.data.frame(sapply(drugTable,prop.table))
# this can't be working the way I want it to. dim(responseProps) is 2 x 90
# I want prop.table of col1 x col2, col1 x col3, col1 x colj, col2 x col3, etc


#testing
test <- as.matrix(drugTable[,1:2])
barplot(test, beside = T, legend = T)
