####################################################################################
####				Tests for SaemixData				####
####################################################################################

# ECO TODO: verifier le nb de sujets ? au moins 2 sujets pour faire tourner saemix ? aussi warning si la reponse n'a que 2 modalites ?

setGeneric(name="read.saemixData",
  def=function(object){standardGeneric("read.saemixData")}
)

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

# SaemixData
source("../R/SaemixData.R")

cat("Beginning tests on SaemixData class...\n")

############ FAIL
##### Missing default values - FAIL with Error
# Creating an empty object: not possible, should print a warning
x<-saemixData()

##### Wrong column names - FAIL
# Wrong Id name : "Creation of saemixData failed"
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="wrongid")

# Wrong predictor name : "Creation of saemixData failed"
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.predictors="wrongpred")

# Wrong response name : "Creation of saemixData failed"
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.response="wrongresp")

# Wrong separator : "Creation of saemixData failed"
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T,sep=",", name.group=1,name.predictors=2,name.response=3)

############ SUCCESS
# No Id name but automatic recognition
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T)

# Column names given as numbers
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group=1,name.predictors=2,name.response=3, name.covariates=4,units=list(x="mg",y="-",covariates="-"))

# Two predictors (silly here)
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.predictors=c("dose","gender"))

# Two predictors, one does not exist - dropped with a warning
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.predictors=c("dose","wrongpred"))

# Wrong number of units for covariates
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(x="mg",y="-",covariates=c("-","-")))
print(x)

# No units for covariates
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(x="mg",y="-"))
print(x)

# No units for x
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(y="-"))
print(x)

# Bimodal response
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.response="gender")

############ SUCCESS
##### Giving all object specifications
x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))
print(x)

# Creating an object, giving only the name - file
x<-saemixData(name.data="../data/PD1.saemix.tab")
print(x)

# Creating an object, giving only the name - data.frame
tab1<-read.table("../data/PD1.saemix.tab",header=T)
x<-saemixData(name.data="tab1")
print(x)

# Creating an object, giving a data.frame 
x<-saemixData(name.data=tab1)

############ SUCCESS
##### Problems with dataset

# X with missing values - remove corresponding lines
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,2]<-dat[41,2]<-dat[55,2]<-NA

x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
print(x)

# Subjects with all missing values - remove entire subject
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[7:9,3]<-NA

x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
print(x)

# Y with missing values - remove corresponding lines
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,3]<-dat[41,3]<-dat[55,3]<-NA

x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
print(x)

# covariate with missing values
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,4]<-dat[41,4]<-dat[55,4]<-NA

x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
summary(x)

# Covariate with character values
dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-vec
x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
print(x)

# Covariate with factor values
dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-as.factor(vec)
x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))
print(x)

############ SUCCESS
##### Method show

x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))
show(x)

showall(x)

##### Method plot

x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))

plot(x)

plot(x,main="Raw data",type="p",col="Blue")

plot(x,main="Raw data",type="b")

plot(x,main="Raw data",type="b",col="Red",ylog=TRUE)

cat("End of tests on SaemixData class.\n")

####################################################################################
####				Tests for SaemixRepData				####
####################################################################################

cat("Beginning tests on SaemixRepData class...\n")

x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))

xrep<-new(Class="SaemixRepData",saemix.data=x)
show(xrep)

cat("End of tests on SaemixRepData class.\n")

####################################################################################
####				Tests for SaemixSimData				####
####################################################################################

cat("Beginning tests on SaemixSimData class...\n")

x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))

xsim<-new(Class="SaemixSimData",saemix.data=x)
show(xsim)

cat("End of tests on SaemixSimData class.\n")
####################################################################################
