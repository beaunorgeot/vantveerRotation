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
# from: https://www.youtube.com/watch?v=V_YNPQoAyCc
# on odds ratio, relative risk, attributal risk (risk difference)
test <- as.matrix(drugTable[,1:2])
barplot(test, beside = T, legend = T) #get 2 plots side by side w/a legend
library(epiR)
#the epi.2by2() returns all 3 measurements if method = cohort.count
epi.2by2(test, method = "cohort.count", conf.level = 0.95)
# Interpreting odds ratio: see 4:30 of video

# The typcial 2x2 format is:row/column, with the exposure (exposuredToDrug1,exposuredToDrug2) as the rows and the outcome as the columns (sens, res)
# yes/yes  yes/no
# no/yes   no/no
# exposure could be gender, or drug. Outcome is some action that exposure can produce. a girl can smoke, a girl can not smoke
# a drug can produce a response/reaction, a drug can produce no response/reaction

# Odds of female smoking are 1.4 times the odds of male smoking
# Odds of being 'yes' to exposure1 are someValue times the odds of 'yes' to exposure2
# If the confidence interval for the odds ratio contains the value of 1, then the results are NOT significant


#I want the joint probability
#For resitance matrix: calculate the number of resistant vs total number in that column.

#make 4 different matrices, 1 for each r to r, r to s, etc.
# Each matrix will be a 90*90. 
rr = matrix(NA, nrow = 90, ncol = 90) #res to drug2 given res to drug1
rs = matrix(NA, nrow = 90, ncol = 90) #res to drug2 given sens to drug1
ss = matrix(NA, nrow = 90, ncol = 90) #sens to drug2 given sens to drug1
sr = matrix(NA, nrow = 90, ncol = 90) #sens to drug2 given res to drug1

for (drug in 1:90){
  for(row in 1:90){
    prior = q3responses[q3responses[,drug]==1,]
    rr[row,drug] = mean(prior[,row])
  }
}

#q3responses %>% filter(Sens_X5.FU == 1) #takes all rows that are sensitive to FU
#ja = q3responses[q3responses[,2]==1,] does the same as the dplyr above
