
setwd("/Users/beaunorgeot/vantveerRotation")
#Headers are in the 3rd row
mainTable <- read.csv("/Users/beaunorgeot/vantveer/SuppTab1.csv", skip = 2)

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
# list to hold dummy vars. Cell lines will either get 0 or 1 for each of these vars
response <- c("Sens","Res") 
# give new column names to the data. rep() for replicate. repeat each item in features x2. repeat each item in response x90. 
# then paste the first item in features to the first item in response using '_' to join them. These + diagnosis become the colnames
drugResponse <- c(paste0(rep(features, 2), "_", rep(response, each=90)))
drugs[,drugResponse] <- NA 


drugs1 <- drugs
# This works!
drugs1 <- drugs1 %>% mutate(Sens_X17.AAG = ifelse(X17.AAG < mean(X17.AAG,na.rm=T),1,0))

# Create Sensitivity DummyVar
for(f in features){
  dummy.var <- paste("Sens", f, sep='_')
  drugs1[, dummy.var] <- NA
  drugs1[,dummy.var] <- ifelse(drugs1[,f] <= mean(drugs1[,f],na.rm=T),1,0) 
}

# Create Resistant DummyVar
for(f in features){
  dummy.var <- paste("Res", f, sep='_')
  drugs1[, dummy.var] <- NA
  drugs1[,dummy.var] <- ifelse(drugs1[,f] >= mean(drugs1[,f],na.rm=T),1,0) 
}


#Test for non-normalcy on the continous outputs
drugsMean <- drugs1 %>% select(1:90)
#acceptable for use w/our data < 2000 samples
#
shtest1 <- lapply(drugsMean,shapiro.test)
#extract statistic and p-values
drugsMeanNormalStats <- sapply(shtest1, `[`, c("statistic","p.value"))
drugsMeanNormalStats <- as.data.frame(drugsMeanNormalStats)
#transpose
drugsMeanNormalStats <- t(drugsMeanNormalStats)
drugsMeanNormalStats <- as.data.frame(drugsMeanNormalStats)
notNormMean <- drugsMeanNormalStats %>% filter(p.value < 0.05) #There are 63 rows that fail the normalcy test
#investigate columns more latter w/qqplots

#Next: sd == sqrt(var) in R, is.normal()?, 3 cuts (low dose is sens, middle 3rd neutral, upper 3rd res)
# check to make sure dummy.var isn't a column, after running correlations on drugs ->
# add the drug target from the other df as a column on this one




#NOTES
# apply(). Write a function that does what you want to 1 column, pass that function and the df of interest to apply()

#1. What is f? mean(drugs1[,f],na.rm=T) creates a logical vector

#EFFECIENT TABLE TRANSFORMATIONS IN R: CONVERT TO MATRIX, DO OPERATIONS, COVERT BACK
#If need store, preallocate as matrix and dump all results into matrix
drugs2 <- as.matrix(drugs1)

#1. column-wise generation of new values
#2. Insertion of new values into matrix, column-wise
#3. convert back to df

# Instead of using mean and sd, perhaps use cut2() to bin each drug into quantiles. So you get
# the lowest 25% become the sensitives, and the highest 25% become the resistants?
#assign cell lines to resistant/sensetive to each agent
#zeroBin <- zero4na %>% mutate_each(ifelse(. < (. - sd(.)),1,0))