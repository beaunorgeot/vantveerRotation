
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
save(inThirds, file = "quartiles3Complete.RData")
q3responses = inThirds[,91:360]
save(q3responses, file = "q3responses.RData") #classification by 3 quantiles
save(drugs1,file = "medUpDown.RData") #classification by median

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
# just create a count summary
sumBidrugs <- bidrugs %>% summarise_each(funs(sum))
sumBidrugs <- as.data.frame(t(sumBidrugs))
sumBidrugs2 <- cbind(ResponseOfDrug = rownames(sumBidrugs), numObservations = sumBidrugs$V1)
#dnames <- names(bidrugs)
#sumBidrugs1 <- cbind(dnames,sumBidrugs)
#colnames(sumBidrugs1) <- c("responseOfDrug", "numObservations ")


#generate correlation table#######
corbidrugs <- cor(bidrugs) #warning: standard deviation is zero
#turn into df so have access to colnames
cors <- as.data.frame(corbidrugs)
inds <- which(abs(cors) > 0.8, arr.ind = T)
length(inds)/2 #there are 182 significant values
dim(cors) #180 by 180. 32400 total pairings. 182/32400 = 0.005617 % of all values that are extremely highly correlated

rnames = rownames(cors)[inds[,1]]
cnames = colnames(cors)[inds[,2]]
# for all columns, if the abs of a column value is >.8, give me the columnName, the rowName, and the actual value ->
# Store these values as a list, that is stored in a larger list

sigValues = vector("list",182) #list() list vs vector in r. Seems NULL generators an atomic vector, which could become a list?

for (i in 1:182){
  sigValues[i] <- c(paste("Value:",sigValues[rnames[i],cnames[i]], "Pairing:", rnames[i], cnames[i], sep=" "))
}

#Error in sigValues[rnames[i], cnames[i]] : incorrect number of dimensions
summary(cors) #there are 5 NAs in each column, maybe this is causing the problems. Could be related to
# the warning that the sd of cors is zero also
##########

#c#### Running Chisq #####
#turning numerics into factors
factoredDrugs <- as.data.frame(lapply(bidrugs,factor))
#summary(factoredDrugs) weired stuff is happening here
freqTable <- sapply(factoredDrugs, table)
countsTable <- as.data.frame(sapply(bidrugs, sum))
# chi1 <- chisq.test(countsTable) not even remotely useful, checks whether all values are associated
###########
#convert to a 2 x j table with sensitivity and resistance as the rows and the drugName as the cols



###### Creating Test case to trouble shoot########
#indexes[,1] is row name
#indexes[,2] is col name
# the following provides the indices of all the values that match. 

bob <- data.frame(replicate(10,sample(-5:5,10,rep=TRUE)))
inds <- which(abs(bob) > 3, arr.ind = T)
#there are 30 postions shown by View(). length() returns 60, a row and a col for each of the 30 values.
# so the number of positions is length(inds)/2. Important for larger sets

#How do I get the values themselves?
rnames = rownames(bob)[inds[,1]]
cnames = colnames(bob)[inds[,2]]

#this seems really close, but isn't working, excpet it seems to do what I want on the last cell
#bob[rnames[1],cnames[1]] :: this works fine
bobsData = vector("list",30)
#bobsData <- NULL
#bobsData = list()
#Below loop does exactly what I want.
for (i in 1:30){
  bobsData[i] <- c(paste("Value:",bob[rnames[i],cnames[i]], "Pairing:",rnames[i], cnames[i], sep=" "))
}
#for (i in cnames) print (cnames[i]) doesn't work
#for (i in 1:30) print (cnames[i]) ::Works.
#for (i in length(cnames)) print (cnames[i]) #prints 1 value

############
test <- bidrugs[,1:5]
chiTest <- sapply(test,chisq.test) #the expected values are all identical and a float
#could just do same thing as chisq.test w/a conditional probability test

##############

#########end test ##################
#Odd's Ratio thoughts for just the first 2 columns of bidrugs. Select just the first 2 columns
# Need row data to create 4 different numbers from those 2 columns

#    n00 = number of cases where x = 0 and y = 0
#    n01 = number of cases where x = 0 and y = 1
#    n10 = number of cases where x = 1 and y = 0
#    n11 = number of cases where x = 1 and y = 1

#the above can be accomplished with table
library(epitools) #for oddsratio()
t <- bidrugs[,1:2]
tp <- table(t)
oddsratio(tp)


# format for Chi square for this experiment
# Row1: Sensitive
# Row 2: Resistant
# Columns 1:90 Each drug