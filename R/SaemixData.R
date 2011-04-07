####################################################################################
####			SaemixData class - definition				####
####################################################################################

###############################
# Definition with initialise
setClass(
  Class="SaemixData",
  representation=representation(
    name.data="character",	# name of dataset
    header="logical",		# for file, whether has header
    sep="character",		# if file, separator
    na="character",		# if file, NA symbol(s)
    name.group="character",	# name of column with ID
    name.predictors="character",# name of column(s) with predictors
    name.response="character",	# name of column with response
    name.covariates="character",# name of column(s) with covariates
    name.X="character",		# name of predictor used on X axis for graphs
    units="list",		# units (list with components for x, y, and cov)
    tab="data.frame",		# data
    N="numeric",		# number of subjects
    index="numeric",		# subject index (id renamed to 1:N)
    id="character",		# subject id
    xind="data.frame",		# matrix of predictors
    y="numeric",		# vector of responses (possibly transformed during fit)
    yorig="numeric",		# vector of responses in original dataset
    cov="data.frame",		# matrix of covariates
    ntot.obs="numeric",		# total number of observations
    nind.obs="numeric"		# number of observations for each subject
  ),
  validity=function(object){
#    cat ("--- Checking SaemixData object ---\n")
    if (length(object@name.data)==0) {
      stop ("[ SaemixData : validation ] Please provide a name for the data (dataset or datafile on disk).")
    }
# Ici ou a la creation, detection automatique ?
    if (nchar(object@name.group)==0) {cat("Missing Id column\n")}
    if (nchar(object@name.predictors[1])==0) {cat("No predictors found\n")}
    if (nchar(object@name.response)==0) {cat("No response found\n")}
    N<-object@N
    if(N!=length(unique(object@id)) | N!=length(unique(object@index)) | N!=length(object@nind.obs)) {
      cat("Size mismatch: id, index and/or nind.obs no longer correspond to the number of subjects.\n")
    }
    if(length(object@tab)>0) {
     if(N<2) {
       cat("Warning: There are only",N,"subjects in the dataset, the SAEM algorithm is a population algorithm designed to analyse longitudinal data from non-linear mixed effect models and may not work with too few subjects.\n")
     }
     if(length(unique(object@y))<3) {
       cat("Warning: The SAEM algorithm currently handles only continuous responses. It seems that the response --",object@name.response,"-- has too few modalities and the statistical model may not be appropriate.\n")
     }
    }
    return(TRUE)
  }
)

setClass(
  Class="SaemixRepData", # Saemix data, replicated for different chains
  representation=representation(
    N="numeric",		# number of subjects
    NM="numeric",		# number of subjects, replicated
    IdM="numeric",		# subject id
    XM="data.frame",		# matrix of predictors
    yM="numeric",		# vector of responses
    nrep="numeric"		# number of replicates
  ),
  validity=function(object){
#    cat ("--- Checking SaemixData object ---\n")
    ntot<-length(object@yM)
    if (length(object@IdM)!=ntot | dim(object@XM)[1]!=ntot) {cat("Size mismatch.\n")}
    return(TRUE)
  }
)

setClass(
  Class="SaemixSimData", # Saemix predicted and simulated data
  representation=representation(
    N="numeric",		# number of subjects
    name.response="character",	# name of column with response
    name.X="character",		# name of predictor used on X axis for graphs
    units="list",		# units (list with components for x, y, and cov)
    index="numeric",		# subject index (id renamed to 1:N)
    id="character",		# subject id
    xind="data.frame",		# matrix of predictors
    y="numeric",		# vector of responses (possibly transformed during fit)
    nsim="numeric",		# number of simulations
    sim.psi="data.frame",	# simulated parameters
    sim.id="numeric",		# id in replications
    sim.rep="numeric",		# replication number
    sim.ypred="numeric",	# simulated predictions (without error)
    sim.y="numeric"		# simulated data (with error)
  ),
  validity=function(object){
#    cat ("--- Checking saemixSimData object ---\n")
    return(TRUE)
  }
)

###############################
# ECO validity ne semble pas etre appele automatiquement quand on cree un objet => il faut l'appeler dans initialize

setMethod(
  f="initialize",
  signature="SaemixData",
  definition= function (.Object,name.data,header,sep,na,name.group, name.predictors,name.response,name.covariates,name.X,units,tab,N,id,y,index, xind,cov,ntot.obs,nind.obs){
#    cat ("--- initialising SaemixData Object --- \n")
    if(missing(name.data)) stop ("Please provide a name for the data (dataset or datafile on disk).")
    .Object@name.data<-name.data
    if(missing(header)) header<-TRUE
    .Object@header<-header
    if(missing(sep)) sep<-""
    .Object@sep<-sep
    if(missing(na)) na<-"NA"
    .Object@na<-na
    if(missing(name.group)) {
      cat("   Missing ID identifier, assuming the ID is in column 1 of the dataset.\n")
      name.group<-"1"
    }
# ECO TODO: reconnaissance automatique (avant affectation a la valeur 2) ?
    if(missing(name.predictors)) {
      name.predictors<-"2"
      cat("   Missing predictors identifier, assuming there is one predictor in column 2 of the dataset.\n")
    }
    if(missing(name.response)) {
      cat("   Missing response identifier, assuming the response is in column 3 of the dataset.\n")
      name.response<-"3"
    }
    if(missing(name.covariates)) name.covariates<-character()
    .Object@name.group<-name.group
    .Object@name.predictors<-name.predictors
    .Object@name.response<-name.response
    .Object@name.covariates<-name.covariates
    if(missing(units)) units<-list(x="-",y="-")
    if(is.null(units$x)) units$x<-"-"
    if(is.null(units$y)) units$y<-"-"
    ncov<-length(name.covariates)
    if(ncov>0) {
      nunit<-length(units$covariates)
      if(nunit==0) units$covariates<-rep("-",ncov)
      if(nunit>ncov) units$covariates<-units$covariates[1:ncov]
      if(nunit<ncov) {
        length(units$covariates)<-ncov
        units$covariates[(nunit+1):ncov]<-"-"
      }
    }
    .Object@units<-units
    .Object@name.X<-name.X
    if(missing(N)) N<-0
    .Object@N<-N
# For completion: these items will be filled by a call to read.saemixData
    if(FALSE) {
      if(missing(tab)) tab<-data.frame()
      .Object@tab<-tab
      .Object@id<-id
      .Object@y<-y
      .Object@index<-index
      .Object@xind<-xind
      .Object@cov<-cov
      .Object@ntot.obs<-ntot.obs
      .Object@nind.obs<-nind.obs
    }
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

# Initialize method for saemixRepData and saemixSimData
setMethod(
  f="initialize",
  signature="SaemixRepData",
  definition= function (.Object,data=NULL,nb.chains=1){
#    cat ("--- initialising SaemixData Object --- \n")
    if(is.null(data)) {
      .Object@N<-.Object@NM<-numeric(0)
      .Object@IdM<-numeric(0)
      .Object@yM<-numeric(0)
      .Object@XM<-data.frame()
    } else {
    N<-data@N
    .Object@N<-N
    .Object@nrep<-nb.chains
    .Object@NM<-N*nb.chains
    IdM<-kronecker(c(0:(nb.chains-1)),rep(N,data@ntot.obs))+rep(data@index,nb.chains)
    yM<-rep(data@y,nb.chains)
    XM<-do.call(rbind,rep(list(data@xind),nb.chains))
    .Object@IdM<-c(IdM)
    .Object@XM<-XM
    .Object@yM<-yM
   }
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

setMethod(
  f="initialize",
  signature="SaemixSimData",
  definition= function (.Object,data=NULL,nsim=1){
#    cat ("--- initialising SaemixData Object --- \n")
    if(!is.null(data)) {
      .Object@N<-data@N
      .Object@id<-data@id
      .Object@index<-data@index
      .Object@name.response<-data@name.response
      .Object@name.X<-data@name.X
      .Object@units<-data@units
      .Object@y<-data@y
      .Object@xind<-data@xind
    }
    .Object@nsim<-nsim
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			SaemixData class - accesseur				####
####################################################################################

# Getteur
setMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.data"={return(x@name.data)},
    "header"={return(x@header)},
    "sep"={return(x@sep)},
    "na"={return(x@na)},
    "name.group"={return(x@name.group)},
    "name.predictors"={return(x@name.predictors)},
    "name.response"={return(x@name.response)},
    "name.covariates"={return(x@name.covariates)},
    "name.X"={return(x@name.X)},
    "units"={return(x@units)},
    "tab"={return(x@tab)},
    "N"={return(x@N)},
    "id"={return(x@id)},
    "y"={return(x@y)},
    "yorig"={return(x@yorig)},
    "index"={return(x@index)},
    "xind"={return(x@xind)},
    "cov"={return(x@cov)},
    "ntot.obs"={return(x@ntot.obs)},
    "nind.obs"={return(x@nind.obs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.data"={x@name.data<-value},
    "header"={x@header<-value},
    "sep"={x@sep<-value},
    "na"={x@na<-value},
    "name.group"={x@name.group<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.response"={x@name.response<-value},
    "name.covariates"={x@name.covariates<-value},
    "name.X"={x@name.X<-value},
    "units"={x@units<-value},
    "tab"={x@tab<-value},
    "N"={x@N<-value},
    "id"={x@id<-value},
    "y"={x@y<-value},
    "yorig"={x@yorig<-value},
    "index"={x@index<-value},
    "xind"={x@xind<-value},
    "cov"={x@cov<-value},
    "ntot.obs"={x@ntot.obs<-value},
    "nind.obs"={x@nind.obs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)


# For saemixRepData
setMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "NM"={return(x@NM)},
    "IdM"={return(x@IdM)},
    "yM"={return(x@yM)},
    "XM"={return(x@XM)},
    "nrep"={return(x@nrep)},
    stop("No such attribute\n")
   )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "NM"={x@NM<-value},
    "IdM"={x@IdM<-value},
    "yM"={x@yM<-value},
    "XM"={x@XM<-value},
    "nrep"={x@nrep<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

# For saemixSimData
setMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "name.response"={return(x@name.response)},
    "name.X"={return(x@name.X)},
    "units"={return(x@units)},
    "id"={return(x@id)},
    "index"={return(x@index)},
    "xind"={return(x@xind)},
    "y"={return(x@y)},
    "nsim"={return(x@nsim)},
    "sim.psi"={return(x@sim.psi)},
    "sim.rep"={return(x@sim.rep)},
    "sim.ypred"={return(x@sim.ypred)},
    "sim.id"={return(x@sim.id)},
    "sim.y"={return(x@sim.y)},
    stop("No such attribute\n")
   )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "name.response"={x@name.response<-value},
    "name.X"={x@name.X<-value},
    "units"={x@units<-value},
    "id"={x@id<-value},
    "index"={x@index<-value},
    "xind"={x@xind<-value},
    "y"={x@y<-value},
    "sim.psi"={x@sim.psi<-value},
    "sim.rep"={x@sim.rep<-value},
    "sim.ypred"={x@sim.ypred<-value},
    "sim.id"={x@sim.id<-value},
    "sim.y"={x@sim.y<-value},
    "nsim"={x@nsim<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
####			SaemixData class - method to read data			####
####################################################################################

setMethod("read.saemixData","SaemixData",
  function(object) {
    ow <- options("warn")
    options("warn"=-1)
# ce test devrait aller dans la definition de la classe
    if(class(object@name.data)!="character") {
    cat("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
    return("Creation of saemixData failed")
  }
    if(exists(object@name.data)) {
      cat("Using the object called",object@name.data,"in this R session as the data.\n")
      dat<-get(object@name.data)
    } else {
      cat("Reading data from file",object@name.data,"\n")
      header<-object@header
      if(is.null(header)) header<-TRUE
      sep<-object@sep
      if(is.null(sep)) sep<-""
      na.strings<-object@na
      if(is.null(na.strings)) na.strings<-"NA"
      dat<-try(read.table(object@name.data,header=header,sep=sep,na.strings=na.strings))
      if(class(dat)=="try-error") stop("The file ",object@name.data," does not exist. Please check the name and path.\n")
      cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
      print(head(dat))
    }
    if(dim(dat)[2]<2) {
      cat("The dataset contains only one column. The non-linear mixed effect model requires at least 3 columns, with subject ID, predictor (at least one) and response. \nPlease check the field separator, currently given as:", paste("sep=\"",object@sep,"\"",sep=""),"\n")
      return("Creation of saemixData failed")
    }
# Automatic recognition of columns
#    ID (one of id, subject or sujet regardless of case)
#    response (one of Y, conc, concentration, resp, response regardless of case)
#    predictors (time and/or dose, regardless of case)
# ECO TODO: improve automatic recognition ?
# check that we have at least a column id, response and X

    if(!is.na(as.integer(object@name.group))) {
# group given as a column number
      object@name.group<-colnames(dat)[as.integer(object@name.group)]
    }
    if(is.na(object@name.group)) object@name.group<-""
    if(object@name.group=="") {
      i1<-match("id",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) {
        i1<-c(grep("subject",tolower(colnames(dat)),fixed=T), grep("sujet",tolower(colnames(dat)),fixed=T))
      }
      if(length(i1)>0) {
        object@name.group<-colnames(dat)[i1[1]]
        cat("    no name for the group variable (ID) given, will use column --",object@name.group,"-- in the dataset.\n")
      }
    }
    if(object@name.group=="" | is.na(match(object@name.group,colnames(dat)))) {
      cat("Please provide a name for the ID column.\n")
      return("Creation of saemixData failed")
    }
   i1<-as.integer(object@name.predictors[!is.na(as.integer(object@name.predictors))])
    if(length(i1)>0) {
      object@name.predictors[!is.na(as.integer(object@name.predictors))]<- colnames(dat)[i1]
    }
    if(is.na(object@name.predictors)) object@name.predictors<-""
    if(length(object@name.predictors)==0 | (length(object@name.predictors)==1 & object@name.predictors[1]=="")) {
      i1<-c(match(c("time","temps","tps","tim","x","dose"),tolower(colnames(dat))))
      i1<-i1[!is.na(i1)]
      if(length(i1)>0) {
        object@name.predictors<-colnames(dat)[i1]
        cat("    no name for the predictor variable given, will use column(s) --",object@name.predictors,"-- in the dataset.\n")
      }
    }
    id1<-match(object@name.predictors,colnames(dat),nomatch=0)
    if(length(id1[id1==0])>0) {
      cat("    cannot find column(s) --",object@name.predictors[id1==0],"-- dropping them from the data.\n")
    }
    xnam<-object@name.predictors[id1>0]
    if(length(xnam)==0) object@name.predictors<-"" else object@name.predictors<-xnam
    if(length(xnam)==0) {
      cat("Please provide at least one predictor.\n")
      return("Creation of saemixData failed")
    }
    if(!is.na(as.integer(object@name.response))) {
# response given as a column number
      object@name.response<-colnames(dat)[as.integer(object@name.response)]
    }
    if(object@name.response=="") {
      i1<-match("y",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) {
        i1<-c(grep("response",tolower(colnames(dat)),fixed=TRUE), match(c("resp","conc"),tolower(colnames(dat))),grep("concentration", tolower(colnames(dat)),fixed=TRUE))
        i1<-i1[!is.na(i1)]
      }
      if(length(i1)>0) {
        object@name.response<-colnames(dat)[i1[1]]
        cat("    no name for the response variable given, will use column --",object@name.response,"-- in the dataset.\n")
      }
    }
    if(is.na(object@name.response)) object@name.response<-""
    if(object@name.response=="" | is.na(match(object@name.response,colnames(dat)))) {
      cat("Please provide a name for the response column.\n")
      return("Creation of saemixData failed")
    }
    if(length(object@name.covariates)>0 & object@name.covariates[1]!="") {
   i1<-as.integer(object@name.covariates[!is.na(as.integer(object@name.covariates))])
      object@name.covariates[!is.na(as.integer(object@name.covariates))]<- colnames(dat)[i1]
    }
    if(nchar(object@name.group)*length(object@name.predictors)* nchar(object@name.response)<=0) {
      stop("Please check the structure of the data file and provide information concerning which columns specify the group structure (ID), the predictors (eg dose, time) and the response (eg Y, conc). See documentation for automatic recognition of column names for these elements.\n")
    }
    if(nchar(object@name.X)==0)
      object@name.X<-object@name.predictors[1]
    if(!is.na(as.integer(object@name.X))) {
      if(dim(dat)[2]<as.integer(object@name.X)) {
        cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
	object@name.X<-object@name.predictors[1]
      } else object@name.X<-colnames(dat)[as.integer(object@name.X)]
    }
    if(match(object@name.X,object@name.predictors,nomatch=0)==0) {
      cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
      object@name.X<-object@name.predictors[1]
    }
# Removing missing values in response & predictor columns
    dat<-dat[!is.na(dat[,object@name.response]),]
    for(i in object@name.predictors) dat<-dat[!is.na(dat[,i]),]
# ECO TODO: what about missing data in covariates & predictor columns
#    for(i in object@name.covariates) dat<-dat[!is.na(dat[,i]),]
    object@tab<-dat

    id<-dat[,object@name.group]
    object@y<-dat[,object@name.response]
    object@id<-as.character(id)
    object@N<-length(unique(object@id))
    xind<-dat[,object@name.predictors,drop=FALSE]
#    if(is.null(dim(xind))) {
#      xind<-data.frame(xind)
#      colnames(xind)<-object@name.predictors
#    }
    object@xind<-xind
    if(length(object@name.covariates)>0) covariates<-dat[,object@name.covariates,drop=FALSE] else covariates<-NULL
#    if(!is.null(covariates) & is.null(dim(covariates))) {
#      covariates<-matrix(covariates,ncol=1)
#      colnames(covariates)<-object@name.covariates
#    }
    object@cov<-as.data.frame(covariates)
    ntot.obs<-length(object@y) # total number of observations
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-nind.obs[match(unique(id),names(nind.obs))]
    object@nind.obs<-c(nind.obs)
    object@index<-rep(1:object@N,times=nind.obs)
    object@ntot.obs<-length(object@y) # total number of observations
#    object@names<-list(group=object@name.group, predictors=object@name.predictors,response=object@name.response, covariates=object@name.covariates)
    options(ow) # reset
    validObject(object)
    return(object)
  }
)

####################################################################################
####			SaemixData class - method to print/show data		####
####################################################################################

setMethod("print","SaemixData",
  function(x,nlines=10,...) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",x@name.data,"\n")
    st1<-paste(x@name.response," ~ ",paste(x@name.predictors,collapse=" + ")," | ", x@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(x@name.predictors)>1) {
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    } else  cat("    Predictor:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    if(length(x@name.covariates)>0) {
      cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
    }
    if(FALSE) {
      cat("    Group column:",x@name.group,"\n")
      cat("    Predictors:",x@name.predictors,"\n")
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
      cat("    Response column:",x@name.response, paste("(",x@units$y,")",sep=""),"\n")
      cat("    Covariates:",x@name.covariates,"\n")
    }
    if(length(x@tab)>0) {
      if(nlines==0) return()
      cat("Dataset characteristics:\n")
      cat("    number of subjects:    ",x@N,"\n")
      if(x@N>0) {
        cat("    number of observations:",x@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(x@nind.obs),digits=digits, nsmall=nsmall), " / ", min(x@nind.obs)," / ",max(x@nind.obs),"\n")
#    if(length(x@tab)>0) print(x@tab)
      }
      if(nlines==(-1)) {
        cat("Data:\n")
        print(x@tab)
      } else {
        cat("First",nlines,"lines of data:\n")
        nrowShow <- min (nlines , nrow(x@tab ))
        print(x@tab[1:nrowShow,])
      }
    } else cat("No data.\n")
  }
)

setMethod("show","SaemixData",
  function(object) {
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(object@name.predictors)>1) {
      cat("    X variable for graphs:",object@name.X, paste("(",object@units$x,")",sep=""),"\n")
    }
    if(length(object@name.covariates)>0) cat("    Covariates:",object@name.covariates,"\n")
    if(length(object@tab)>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@tab ))
      print(object@tab[1:nrowShow,])
    } else cat("No data.\n")
  }
)

# Could be print, with only head of data
setMethod("showall","SaemixData",
  function(object) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    cat("    header:",object@header,"\n")
    cat("    sep:",object@sep,"\n")
    cat("    na:",object@na,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    cat("    subject identifier:    ",object@name.group,"\n")
    cat("    predictors:       ",object@name.predictors,"\n")
    cat("    response:         ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
    cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
    if(length(object@name.covariates)>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
    }
    cat("Dataset characteristics:\n")
    cat("    number of subjects:    ",object@N,"\n")
    if(object@N>0) {
      cat("    number of observations:",object@ntot.obs,"\n")
      cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@tab)>0) print(object@tab)
    }
    if(length(object@tab)>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@tab ))
      ncolShow <- min (10 , ncol(object@tab))
      print(object@tab[1:nrowShow,])
    } else cat("No data.\n")
  }
)

# SaemixRepData
setMethod("show","SaemixRepData",
  function(object) {
    cat("Object of class saemixRepData\n")
    cat("    replicated data used in the SAEM algorithm\n")
    cat("    number of subjects in initial dataset",object@N,"\n")
    cat("    number of replications",object@nrep,"\n")
    cat("    number of subjects in replicated dataset",object@NM,"\n")
  }
)

# SaemixSimData
setMethod("show","SaemixSimData",
  function(object) {
    cat("Object of class SaemixSimData\n")
    cat("    data simulated according to a non-linear mixed effect model\n")
    cat("Characteristics of original data\n")
    cat("    number of subjects:",object@N,"\n")
    cat("    summary of response:\n")
    print(summary(object@y))
    cat("Characteristics of simulated data\n")
    if(length(object@sim.y)>0) {
      cat("    number of simulated datasets:",object@nsim,"\n")
      cat("    summary of simulated response\n")
      print(summary(object@sim.y))
    } else cat("    no simulations performed yet\n")
# ECO TODO: decide if predictions are in results or in a data-type object
    if(FALSE) {
    cat("Characteristics of model predictions\n")
    if(length(object@ppred)>0) {
      cat("    summary of population predictions\n")
      print(summary(object@ppred))
    } else cat("    population predictions not obtained yet\n")
    if(length(object@ypred)>0) {
      cat("    summary of mean predictions\n")
      print(summary(object@ypred)>0)
    } else cat("    mean predictions not obtained yet\n")
    if(length(object@ipred)>0) {
      cat("    summary of individual predictions\n")
      print(summary(object@ipred))
    } else cat("    individual predictions not obtained yet\n")
    }
# ECO TODO; autres resumes, y compris tests sur les npde, moyennes npde, pd, etc...
  }
)

####################################################################################
####				Summary method for SaemixData			####
####################################################################################

setMethod("summary","SaemixData",
  function(object) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(object@name.predictors)>1) cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
    if(length(object@name.covariates)>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
    }
    cat("Dataset characteristics:\n")
    cat("    number of subjects:    ",object@N,"\n")
    if(object@N>0) {
      cat("    number of observations:",object@ntot.obs,"\n")
      cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@tab)>0) print(object@tab)
    }
    res<-list(N=object@N,nobs=list(ntot=object@ntot.obs,nind=object@nind.obs), id=object@id,xind=object@xind,covariates=object@cov,y=object@y)
    invisible(res)
 }
)

####################################################################################
####			SaemixData class - method to plot			####
####################################################################################

saemix.data.setoptions<-function(saemix.data) {
# setting default plot options
  plot.opt<-list(
# Options for plot types
    ilist=c(1:saemix.data["N"]),
    level=0:1,
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main="",				# title
    xlab="",
    ylab="",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1)

     if(is.null(plot.opt$name.X))
        plot.opt$name.X<-saemix.data["name.predictors"][1]
    plot.opt$xlab<-paste(plot.opt$name.X," (",plot.opt$units$x,")", sep="")
     if(length(saemix.data["name.response"])>0)
    plot.opt$ylab<-paste(saemix.data["name.response"]," (",saemix.data["units"]$y,")", sep="")
   return(plot.opt)
}

replace.data.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  if(length(args1)>2) {
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0)
    plot.opt[[names(args1)[i]]]<-eval(args1[[i]]) else {
      if(names(args1)[i]!="plot.type") cat("Argument",names(args1)[i],"not available, check spelling.\n")
    }
   }
  }
  return(plot.opt)
}

# Plot the data, either as points or as lines grouped by x@name.group
setMethod("plot","SaemixData",
  function(x,y,...) {
    if(length(x@tab)>0) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"l"
    plot.opt<-saemix.data.setoptions(x)
    plot.opt$new<-TRUE
    plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
    plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
    plot.opt<-replace.data.options(plot.opt,...)
    logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
    if(plot.opt$new) par(mfrow=c(1,1))
      if(plot.type=="p" | plot.type=="b") {
        plot(x@tab[,x@name.X],x@tab[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
      if(plot.type=="l") {
        plot(x@tab[,x@name.X],x@tab[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
      }
      if(plot.type=="l" | plot.type=="b") {
        for(isuj in unique(x@tab[,x@name.group])) {
          lines(x@tab[x@tab[,x@name.group]==isuj,x@name.X], x@tab[x@tab[,x@name.group]==isuj,x@name.response],col=plot.opt$col, lty=plot.opt$lty,lwd=plot.opt$lwd)
      }
      }
    } else cat("No data to plot.\n")
  }
)

# Check for mirror plots
setMethod("plot","SaemixSimData",
  function(x,y,irep=-1,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"l"
    plot.opt<-saemix.data.setoptions(x)
    plot.opt$new<-TRUE
    plot.opt$plot.type<-"b"
    plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
    plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
    plot.opt<-replace.data.options(plot.opt,...)
    logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
    if(length(x@sim.y)==0) cat("No simulated data.\n") else {
      if(irep<0) irep<-sample(unique(x@sim.rep),1)
      tit<-paste("Mirror plot (replication ",irep,")",sep="")
      tab<-cbind(id=x@id,x=x@xind[,x@name.X],y=x@sim.y[x@sim.rep==irep])
      if(plot.type=="p" | plot.type=="b") {
        plot(tab[,"x"],tab[,"y"],xlab=plot.opt$xlab, ylab=plot.opt$ylab, col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=tit,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
      if(plot.type=="l") {
        plot(tab[,"x"],tab[,"y"],type="n",xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=tit, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
      }
      if(plot.type=="l" | plot.type=="b") {
        for(isuj in unique(tab[,"id"])) {
          lines(tab[tab[,"id"]==isuj,"x"],tab[tab[,"id"]==isuj,"y"])
      }
      }
    }
  }
)

####################################################################################
####			SaemixData class - User-level function			####
####################################################################################

saemixData<-function(name.data,header,sep,na,name.group, name.predictors,name.response,name.X,name.covariates=c(), units=list(x="",y="",covariates=c())) {
# setting proper types for the SaemixData class
  if(missing(name.data)) {
    cat("Error in saemixData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of SaemixData failed")
  }
  if(is.data.frame(name.data)) name.data<-deparse(substitute(name.data))
  if(missing(header)) header<-TRUE
  if(missing(sep)) sep<-""
  if(missing(na)) na<-"NA" else {na<-as.character(na);na[is.na(na)]<-"NA"}
  if(missing(name.group)) name.group<-"" else name.group<-as.character(name.group)
  if(missing(name.predictors)) name.predictors<-"" else name.predictors<-as.character(name.predictors)
  if(missing(name.response)) name.response<-"" else  name.response<-as.character(name.response)
  if(missing(name.X)) name.X<-"" else name.X<-as.character(name.X)
  name.covariates<-as.character(name.covariates)
  x<-new(Class="SaemixData",name.data=name.data,header=header,sep=sep,na=na, name.group=name.group,name.predictors=name.predictors,name.response=name.response, name.covariates=name.covariates,units=units,name.X=name.X)
#  showall(x)
  x1<-read.saemixData(x)
  if(class(x1)!="character") cat("\n\nThe following SaemixData object was successfully created:\n\n")
  print(x1,nlines=0)
  return(x1)
}

####################################################################################
