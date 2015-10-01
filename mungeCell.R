
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


#Test for non-normalcy on the continous outputs
drugsMean <- drugs1 %>% select(1:90)

shtest1 <- lapply(drugsMean,shapiro.test)
#extract statistic and p-values
drugsMeanNormalStats <- sapply(shtest1, `[`, c("statistic","p.value"))
drugsMeanNormalStats <- as.data.frame(drugsMeanNormalStats)
#transpose
drugsMeanNormalStats <- t(drugsMeanNormalStats)
drugsMeanNormalStats <- as.data.frame(drugsMeanNormalStats)
drugsMeanNormalStats <- add_rownames(drugsMeanNormalStats,"Drug Name")
#acceptable for use w/our data < 2000 samples
# If the chosen alpha level is 0.05 and the p-value is less than 0.05, 
#then the null hypothesis that the data are normally distributed is rejected
notNormMean <- drugsMeanNormalStats %>% filter(p.value < 0.05) #There are 63 rows that fail the normalcy test
#This might mean that using anything related to sd isn't a good idea
# Consider just using quantiles
#investigate columns more latter w/qqplots

#Next: sd == sqrt(var) in R, 3 cuts (low dose is sens, middle 3rd neutral, upper 3rd res)
#QUANTILES
drugsCut <- drugs
drugsCutNames <- names(drugsCut)

#impute w/median: try on single column
#works fine
#drugsCut$X17.AAG <- ifelse(is.na(drugsCut$X17.AAG), median(drugsCut$X17.AAG, na.rm= TRUE), drugsCut$X17.AAG)
#colnames(drugsCut)[1] #colname here needs the quotes
#drugsCut[,"X17.AAG"] <- ifelse(is.na(drugsCut[,"X17.AAG"]), median(drugsCut[,"X17.AAG"], na.rm= TRUE), drugsCut[,"X17.AAG"])

#programatic works
#for (f in drugsCutNames) {
#  drugsCut[,f] <- ifelse(is.na(drugsCut[,f]), median(drugsCut[,f], na.rm= TRUE), drugsCut[,f])
#}

#Doesn't work: check me. I think it just needs drugsCut$X17.AAG instead of X17.AAG. Remember that is.na() is interpreting it, not dplyr
#drugsCut <- drugsCut %>% mutate(ifelse(is.na(X17.AAG),0,X17.AAG))

library(Hmisc)
inThirds <- as.data.frame(lapply(drugsCut, cut2, g=3))
thirdNames <- names(inThirds)
#this produces a df that looks correct
for (n in thirdNames) {levels(inThirds[,n]) <- c('Sens','AVG','Res')}
# but: levels(inThirds) returns NULL
#this fixes them
levels(inThirds) <- c('Sens','AVG','Res')
#test[,"X17.AAG"] <- ifelse(test[,"X17.AAG"] == "Sens",1,0) WORKS!

#Sens dummy
for (n in thirdNames) {
  dummy.var <- paste("Sens", n, sep='_')
  inThirds[,dummy.var] <- NA
  inThirds[,dummy.var] <- ifelse(inThirds[,n] == "Sens",1,0)
}

# Res dummy
for (n in thirdNames) {
  dummy.var <- paste("Res", n, sep='_')
  inThirds[,dummy.var] <- NA
  inThirds[,dummy.var] <- ifelse(inThirds[,n] == "Res",1,0)
}

# Avg dummy
for (n in thirdNames) {
  dummy.var <- paste("AVG", n, sep='_')
  inThirds[,dummy.var] <- NA
  inThirds[,dummy.var] <- ifelse(inThirds[,n] == "AVG",1,0)
}

# BRING IN DRUG TARGET DATA
targIn <- read.csv("~/vantveerRotation/targetClinicalModel.csv",  skip = 2)
# just get drug name and target name, and change the names
targOnly <- targIn %>% select(drugName = Drug.compound,target = Target)
targOnly <- targOnly[1:93,] #no data below r93
targOnly <- targOnly[-c(53,54),] #empty rows
targOnly$drugName <- as.character(targOnly$drugName)

#NEXT: 
#1.Remove drugs that aren't in the drugs df, and make sure all that are in the drugs df are present in targOnly$drugName
#2. Create X17.AAG_target column and add in the target?

#add cell descriptions back in
allThirds <- cbind(cellDesciptions,inThirds)
alldrugs <- cbind(cellDesciptions,drugs1)

#save working dfs to data file
save(inThirds, file = "quartiles3.RData")
save(drugs1,file = "medUpDown.RData")

#Next: 
# after running correlations on drugs ->add the drug target from the other df as a column on this one

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

# this all belongs in exploreCell
#CORRELATION
#take only the dummy var responses so that we're in the world of Sens and Res
bidrugs <- drugs1[,91:270]
#generate correlation table
corbidrugs <- cor(bidrugs) #warning: standard deviation is zero
#turn into df so have access to colnames
c <- as.data.frame(corbidrugs)
#What correlation coeffecient values are highly signficant?