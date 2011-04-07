####################################################################################
####		Defining class containing data, model and options		####
####################################################################################

###############################
# definition
setClass(Class="SaemixObject",
  representation=representation(
    data="SaemixData",		# Data
    model="SaemixModel",	# Model
    results="SaemixRes",	# Fit results
    rep.data="SaemixRepData",	# Data replicates during algorithm (nb chains)
    sim.data="SaemixSimData", 	# Simulated data
    options="list",		# Options and parameters for algorithm
    prefs="list"		# Options for graphs
  ),
  validity=function(object){
#    cat ("--- Checking SaemixObject object ---\n")
    validObject(object@data)
    validObject(object@model)
    return(TRUE)
  }
)

###############################
# Initialize
setMethod(
  f="initialize",
  signature="SaemixObject",
  definition= function (.Object,data,model,options=list()){
#    cat ("--- initialising SaemixObject --- \n")
    .Object@data<-data
    .Object@model<-model
# ECO TODO: tidy up (check omega.init defined in call to Model)
if(FALSE){
    omega.init<-model@omega.init
    indest.omega<-which(model@covariance.model>0)
    i1.omega2<-which(diag(model@covariance.model)>0)
    i0.omega2<-which((1-diag(model@covariance.model))>0)
    omega.init[i0.omega2,i0.omega2]<-diag(length(i0.omega2))
    .Object@model@omega.init<-omega.init
    .Object@model@i0.omega2<-i0.omega2
    .Object@model@i1.omega2<-i1.omega2
    .Object@model@indest.omega<-indest.omega
    if(dim(model@covariate.model)[1]>length(data@name.covariates)) {
      cat("Number of covariates in model (",dim(model@covariate.model)[1],") is larger than number of covariates in the dataset (",length(data@name.covariates),"), keeping only the first.",length(data@name.covariates),"\n")
      model@covariate.model<-model@covariate.model[1:length(data@name.covariates),]
    }
    if(dim(model@covariate.model)[1]<length(data@name.covariates)) {
      stop("The number of lines in the covariate model needs to be equal to the number of covariates\n")
    }
    covariate.model<-matrix(c(rep(1,model@nb.parameters), c(t(model@covariate.model))),ncol=model@nb.parameters,byrow=TRUE)
# setting the names of the fixed effects
    if(dim(model@covariate.model)[1]>0) {
      rownames(covariate.model)<-c("Fixed",data@name.covariates)
      nam.with.cov<-rep("",length(model@covariate.model))
      row1<-matrix(rep(model@name.modpar,length(data@name.covariates)), ncol=length(model@name.modpar),byrow=TRUE)
      col1<-matrix(rep(data@name.covariates,length(model@name.modpar)), ncol=length(model@name.modpar))
      idcov<-which(model@covariate.model==1)
      nam.with.cov[idcov]<-paste("beta_",col1[idcov],"(",row1[idcov],")",sep="")
      nam1<-rbind(model@name.modpar,matrix(nam.with.cov, ncol=length(model@name.modpar)))
      nam1<-c(nam1)
      .Object@model@name.fixed<-nam1[nam1!=""]
    } else {
      rownames(covariate.model)<-c("Fixed")
      .Object@model@name.fixed<-model@name.modpar
    }
#covariate.model<-rbind(rep(1,nb.parameters),model@covariate.model)
    .Object@model@covariate.model<-covariate.model
}
    .Object@model@name.cov<-.Object@data@name.covariates
    if(length(.Object@model@name.cov)>0 & sum(.Object@model@covariate.model)>0) {
      try(rownames(.Object@model@covariate.model)<-.Object@model@name.cov)
      try(rownames(.Object@model@betaest.model)[2:(1+ length(.Object@model@name.cov))]<-.Object@model@name.cov)
    }
    opt<-saemixControl()
    if(length(options)>0) {
    for(i in names(options)) opt[i]<-options[i]
      if(!opt$fix.seed) {
      rm(.Random.seed)
      runif(1)
      opt$seed<-.Random.seed[5]
      }
    if(is.null(options$nb.chains) & data@N>0) opt$nb.chains<-ceiling(50/data@N)
    if(data@N>0 && data@N<50 & opt$nb.chains<ceiling(50/data@N)) {
      cat("The number of subjects is small, increasing the number of chains to", ceiling(50/data@N),"to improve convergence\n")
      opt$nb.chains<-ceiling(50/data@N)
    }
    if(opt$ipar.lmcmc<2) {
      opt$ipar.lmcmc<-2
      cat("Value of L_MCMC too small, setting it to 2 (computation of the conditional means and variances of the individual parameters)\n")
    }
    }
    opt$nbiter.sa<-round(opt$nbiter.saemix[1]/2)
    opt$nbiter.tot<-sum(opt$nbiter.saemix)
    .Object@options<-opt
    .Object@prefs<-saemix.plot.setoptions(.Object)
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			saemixObject class - accesseur				####
####################################################################################

# Getteur
setMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "data"={return(x@data)},
    "model"={return(x@model)},
    "results"={return(x@results)},
    "rep.data"={return(x@rep.data)},
    "sim.data"={return(x@sim.data)},
    "options"={return(x@options)},
    "prefs"={return(x@prefs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "data"={x@data<-value},
    "model"={x@model<-value},
    "results"={x@results<-value},
    "rep.data"={x@rep.data<-value},
    "sim.data"={x@sim.data<-value},
    "options"={x@options<-value},
    "prefs"={x@prefs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)


####################################################################################
####				Summary method for SaemixObject			####
####################################################################################
setMethod("summary","SaemixObject",
  function(object) {
    if(length(object@results@fixed.effects)==0) {
      cat("Object of class SaemixObject, no fit performed yet.\n")
      return()
    }
    digits<-2;nsmall<-2
    cat("----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@results@se.fixed)==0) {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.res[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]))
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.res[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]),c(object@results@se.fixed,object@results@se.respar[object@results@indx.res]), stringsAsFactors=FALSE)
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(object@results@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@results@indx.cov]<-1-normcdf(abs(wstat[object@results@indx.cov]))
      tab<-cbind(tab,"p-value"=pval,stringsAsFactors=FALSE)
      }
    }
    tab.fix<-tab
    for(i in 2:dim(tab)[2]) {
     xcol<-as.double(as.character(tab[,i]))
     idx<-which(!is.na(xcol) & xcol!="-")
     tab[idx,i]<-format(xcol[idx],digits=digits,nsmall=nsmall)
    }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    if(length(object@results@se.omega2)==0) {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega], object@results@se.omega2[object@results@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
    tab.random<-tab
    for(i in 2:dim(tab)[2])
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("------  Correlation matrix of random effects  ------\n")
    cat("----------------------------------------------------\n")
    tab<-cov2cor(object@results@omega[object@results@indx.omega, object@results@indx.omega])
    tab.corr<-tab
    for(i in 1:dim(tab)[2])
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    try(colnames(tab)<-rownames(tab)<-object@results@name.random)
    print(tab,quote=FALSE)
    l1<-rep(NA,3)
    tab.ll<-data.frame(Method=c("Linearisation","Importance Sampling","Gaussian Quadrature"),"-2xLL"=l1,AIC=l1,BIC=l1)
    if(length(object@results@ll.lin)>0 | length(object@results@ll.is)>0 | length(object@results@ll.gq)>0) {
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@results@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@results@ll.lin),"\n")
    cat("      AIC =",object@results@aic.lin,"\n")
    cat("      BIC =",object@results@bic.lin,"\n")
    tab.ll[1,2:4]<-c((-2*object@results@ll.lin),object@results@aic.lin, object@results@bic.lin)
    }
    if(length(object@results@ll.is)>0) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@results@ll.is),"\n")
    cat("      AIC =",object@results@aic.is,"\n")
    cat("      BIC =",object@results@bic.is,"\n")
    tab.ll[2,2:4]<-c((-2*object@results@ll.is),object@results@aic.is, object@results@bic.is)
    }
    if(length(object@results@ll.gq)>0) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@results@ll.gq),"\n")
    cat("      AIC =",object@results@aic.gq,"\n")
    cat("      BIC =",object@results@bic.gq,"\n")
    tab.ll[3,2:4]<-c((-2*object@results@ll.gq),object@results@aic.gq, object@results@bic.gq)
    }
    cat("----------------------------------------------------\n")
    }
    tab<-data.frame(Id=unique(object@data@id),object@results@cond.mean.psi)
    try(colnames(tab)[-c(1)]<-object@model@name.modpar,silent=TRUE)
    npar<-length(object@results@name.fixed)
    coef<-list(fixed=tab.fix[1:npar,2],random=list(map.psi=object@results@map.psi, cond.mean.psi=tab))
    sigma<-tab.fix[-c(1:npar),2]

    res<-list(fixed.effects=tab.fix,sigma=sigma,random.effects=tab.random, correlation.matrix=tab.corr,logLik=tab.ll,coefficients=coef)
    if(length(object@results@fim)>0) res$FIM<-object@results@fim
    res$data<-list(N=object@data@N,nobs=list(ntot=saemix.fit@data@ntot.obs, nind=saemix.fit@data@nind.obs),id=object@data@id,xind=object@data@xind, y=object@data@y)
    if(length(object@results@ypred)>0 | length(object@results@ipred)>0  | length(object@results@ppred)>0 | length(object@results@icpred)>0) {
      res$fitted<-list(population=list(pop.param=object@results@ppred, pop.mean=object@results@ypred),individual=list(map.ipred=object@results@ipred, cond.ipred=object@results@icpred))
    }
     if(length(object@results@wres)>0 | length(object@results@iwres)>0  | length(object@results@icwres)>0 | length(object@results@pd)>0) {
      res$residuals<-list(population=list(wres=object@results@wres), individual=list(map.iwres=object@results@iwres,cond.iwres=object@results@icwres, pd=object@results@pd, npde=object@results@npde))
    }

    invisible(res)
 }
)

####################################################################################
####			Print and show methods for SaemixObject			####
####################################################################################

setMethod("print","SaemixObject",
  function(x,nlines=10,...) {
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    print(x@data,nlines=nlines)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    print(x@model)
    cat("-----------------------------------\n")
    cat("----    Key algorithm options  ----\n")
    cat("-----------------------------------\n")
    algs<-c("MAP","FIM","LL by IS")
    if(sum(x@options$algorithms)>0)
      cat("    Algorithms:",paste(algs[x@options$algorithms==1],collapse=", "),"\n") else cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),x@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",x@options$nb.chains,"\n")
    cat("    Seed: ",x@options$seed,"\n")
    cat("    Number of MCMC iterations for IS: ",x@options$nmc.is,"\n")
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",x@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",x@options$nb.simpred,"\n")
    cat("    Input/output\n")
    cat("        save the results to a file: ",x@options$save,"\n")
    cat("        save the graphs to files: ",x@options$save.graphs,"\n")
    cat("        directory where results should be saved: ",x@options$directory,"\n")
    cat("-----------------------------------\n")
    cat("----          Results          ----\n")
    cat("-----------------------------------\n")
    print(x@results)
  }
)

setMethod("show","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------------\n")
    cat("----         Data and Model          ----\n")
    cat("-----------------------------------------\n")
#    show(object@data)
    cat("Data\n")
    cat("    Dataset",object@data@name.data,"\n")
    st1<-paste(object@data@name.response," ~ ",paste(object@data@name.predictors, collapse=" + ")," | ", object@data@name.group,sep="")
    cat("    Longitudinal data:",st1,"\n\n")
#    show(object@model)
    cat("Model:\n")
    if(length(object@model@description)>0 && nchar(object@model@description)>0) {
      cat("   ",object@model@description,"\n")}
    fix1<-ifelse(object@model@fixed.estim==1,""," [fixed]")
    cat("    ",object@model@nb.parameters,"parameters:", paste(object@model@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@model@error.model,"\n")
    if(dim(object@model@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@model@covariate.model)
    } else cat("     No covariate\n")
    cat("\n")
    cat("Key options\n")
    algs<-c("MAP","FIM","LL by IS")
    if(sum(object@options$algorithms)>0)
      cat("    Algorithms:",paste(algs[object@options$algorithms==1],collapse=", "),"\n") else cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Seed: ",object@options$seed,"\n")
    cat("    Number of MCMC iterations for IS: ",object@options$nmc.is,"\n")
    cat("    Input/output\n")
    if(object@options$save)
      cat("        save the results to a file: ",object@options$save,"\n") else
      cat("        no graphs\n")
    if(object@options$save.graphs)
      cat("        save the graphs to files: ",object@options$save.graphs,"\n") else
      cat("        no graphs\n")
    if(object@options$save | object@options$save.graphs)
      cat("        directory where results are saved: ",object@options$directory,"\n")
    if(FALSE) {
    if(length(object@rep.data)>0) irep<-1 else irep<-0
    if(length(object@sim.data)>0) isim<-1 else isim<-0
    if(irep>0 | isim>0) {
      cat("-----------------------------------\n")
      cat("----      Other components     ----\n")
      cat("-----------------------------------\n")
      if(irep>0) cat("    Replicated data on ",object@options$nb.chains," chains\n")
      if(isim>0) cat("    Simulated data, ",object@options$nb.sim," simulations\n")
    }
    cat("-----------------------------------\n")
  }
  if(length(object@results@fixed.effects)>0) {
    cat("-----------------------------------------\n")
    cat("----              Results            ----\n")
    cat("-----------------------------------------\n")
    show(object@results)
  }}
)

# Could be print, with only head of data
setMethod("showall","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    showall(object@data)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    show(object@model)
    cat("-----------------------------------\n")
    cat("----      Algorithm options    ----\n")
    cat("-----------------------------------\n")
    algs<-c("MAP","FIM","LL by IS")
    if(sum(object@options$algorithms)>0)
      cat("    Algorithms:",paste(algs[object@options$algorithms==1],collapse=", "),"\n") else cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Number of iterations: ",st1,"\n")
    cat("        nb of iterations for SA: ",object@options$nbiter.sa,"\n")
    cat("        nb of burning iterations: ",object@options$nbiter.burn,"\n")
    cat("    Seed:\n")
    cat("        setting a random seed: ",!object@options$fix.seed,"\n")
    cat("        seed for the random number generator: ",object@options$seed,"\n")
    cat("    Estimation of LL by Importance Sampling:\n")
    cat("        number of MCMC samples: ",object@options$nmc.is,"\n")
    cat("        nu for IS: ",object@options$nu.is,"\n")
    cat("        produce plots during the estimation of LL by IS: ",object@options$print.is,"\n")
    cat("    Input/output\n")
    cat("        display progress during the estimation process: ",object@options$displayProgress,"\n")
    cat("        nb of iterations after which to display progress: ",object@options$nbdisplay,"\n")
    cat("        print out the results after the fit: ",object@options$print,"\n")
    cat("        save the results to a file: ",object@options$save,"\n")
    cat("        save the graphs to files: ",object@options$save.graphs,"\n")
    cat("        directory where results should be saved: ",object@options$directory,"\n")
    cat("        whether warnings should be output during the fit: ",object@options$warnings,"\n")
    cat("    SAEM algorithm\n")
    cat("        number of MCMC iterations for each kernel: ",object@options$nbiter.mcmc,"\n")
    cat("        probability of acceptance: ",object@options$proba.mcmc,"\n")
    cat("        : ",object@options$stepsize.rw,"\n")
    cat("        : ",object@options$rw.init,"\n")
    cat("        : ",object@options$alpha.sa,"\n")
    cat("        maximum nb of iterations for estimation of fixed effects: ",object@options$maxim.maxiter,"\n")
    cat("    Estimation of LL by Gaussian Quadrature:\n")
    cat("        number of nodes: ",object@options$nnodes.gq,"\n")
    cat("        width of integral: ",object@options$nsd.gq,"\n")
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",object@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",object@options$nb.simpred,"\n")
    cat("    Estimation of individual parameters\n")
    cat("        nb of iterations: ",object@options$ipar.lmcmc,"\n")
    cat("        : ",object@options$ipar.rhomcmc,"\n")
    cat("        size of confidence interval: ",object@options$ipar.rmcmc,"\n")
    cat("-----------------------------------\n")
    if(length(object@results@fixed.effects)>0) {
      cat("-----------------------------------------\n")
      cat("----              Results            ----\n")
      cat("-----------------------------------------\n")
      show(object@results)
    }
  }
)

####################################################################################
####			SaemixObject class - method to predict			####
####################################################################################

setMethod(f="predict",
  signature="SaemixObject",
  def=function(object,newdata,...) {
    if(length(object["results"]["map.psi"])==0)
      object<-map.saemix(object)
# en principe n'arrive jamais car on les calcule pdt le fit...
    if(length(object["results"]["cond.mean.phi"])==0)
      object<-conddist.saemix(object)
    saemix.res<-object["results"]
# Individual predictions
    ipred<-object["model"]["model"](saemix.res["map.psi"][, 2:dim(saemix.res["map.psi"])[2]],object["data"]["index"], object["data"]["xind"])
    psiM<-transphi(saemix.res["cond.mean.phi"],object["model"]["transform.par"])
    icond.pred<-object["model"]["model"](psiM,object["data"]["index"], object["data"]["xind"])
    saemix.res["ipred"]<-ipred
    saemix.res["icpred"]<-icond.pred
# Individual weighted residuals
    ares<-saemix.res["respar"][1]
    bres<-saemix.res["respar"][2]
    gpred<-cutoff(ares+bres*abs(ipred))
    iwres<-(ipred-object["data"]["y"])/gpred
    gpred<-cutoff(ares+bres*abs(icond.pred))
    icwres<-(icond.pred-object["data"]["y"])/gpred
    saemix.res["iwres"]<-iwres
    saemix.res["icwres"]<-icwres
# Population predictions using the population parameters [ f(mu) ]
    psiM<-transphi(saemix.res["mean.phi"],object["model"]["transform.par"])
    ppred<-object["model"]["model"](psiM,object["data"]["index"], object["data"]["xind"])
    saemix.res["ppred"]<-ppred
# Population weighted residuals: needs the individual variance-covariance matrix => use compute.sres to estimate these by simulations
    object["results"]<-saemix.res
    return(object)
  }
)

####################################################################################
####			SaemixObject class - method to plot			####
####################################################################################

setMethod(f="plot",
  signature="SaemixObject",
  def=function(x,y,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("plot.type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"reduced"
#    cat("plot.type=",plot.type,"\n")
    if(plot.type[1]=="reduced") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions")
    if(plot.type[1]=="full") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions","residuals.scatter","residuals.distribution","vpc")

    pltyp<-c("data","convergence","likelihood","individual.fit", "population.fit", "both.fit","observations.vs.predictions","residuals.scatter", "residuals.distribution","vpc","npde","random.effects","marginal.distribution", "correlations","parameters.vs.covariates","randeff.vs.covariates")
    ifnd<-pmatch(plot.type,pltyp)
    if(sum(is.na(ifnd))>0) {
      cat("The following plot types were not found or are ambiguous:", plot.type[is.na(ifnd)],"\n")
    }
    ifnd<-ifnd[!is.na(ifnd)]
    if(length(ifnd)==0) return()
    plot.type<-pltyp[ifnd]
    interactive<-x["prefs"]$interactive
    id.pred<-match(plot.type,c("observations.vs.predictions","individual.fit", "both.fit","residuals.scatter","residuals.distribution"))
    id.map<-match(plot.type,c("randeff.vs.covariates"))
    id.sim<-match(plot.type, c("vpc","npde","residuals.scatter","residuals.distribution"))
    id.pred<-id.pred[!is.na(id.pred)]
    id.sim<-id.sim[!is.na(id.sim)]
    id.map<-id.map[!is.na(id.map)]
    namObj<-deparse(substitute(x))
#    cat(namObj,"\n")
    if(length(id.pred)>0) {
      if(length(x["results"]["ipred"])==0 | length(x["results"]["iwres"])) {
        boolpred<-TRUE
        if(interactive) {
      cok<-readline(prompt="Computations will be performed to obtain model predictions, proceed ? (y/Y) [default=yes] ")
        if(cok!="y"&cok!="Y"&cok!="yes"&cok!="") boolpred<-FALSE}
        if(boolpred) {
          x<-predict(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.sim)>0) {
      if(length(x["results"]["npde"])==0) {
        boolpred<-TRUE
        if(interactive) {
        cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
        if(cok!="y"&cok!="Y"&cok!="yes"&cok!="") boolpred<-FALSE}
        if(boolpred) {
          x<-compute.sres(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.map)>0) {
      if(length(x["results"]["map.eta"])==0) {
        cat("Computing ETA estimates and adding them to fitted object.\n")
	x<-compute.eta.map(x)
        assign(namObj,x,envir=parent.frame())
      }
    }
    for(ipl in plot.type) {
      switch (EXPR=ipl,
    "data"={
       cat("Plotting the data\n")
       saemix.plot.data(x,...)
    },
    "convergence"={
       cat("Plotting convergence plots\n")
       saemix.plot.convergence(x,...)
    },
    "likelihood"={
       cat("Plotting the likelihood\n")
       saemix.plot.llis(x,...)
    },
    "observations.vs.predictions"={
      if(length(x["results"]["ipred"])>0) {
        cat("Plotting observations versus predictions\n")
        saemix.plot.obsvspred(x,...)
      }
    },
    "individual.fit"={
      x@prefs$level<-c(1)
      if(length(x["results"]["ipred"])>0) {
        cat("Plotting individual fits\n")
        saemix.plot.fits(x,...)
      }
    },
    "population.fit"={
      x@prefs$level<-c(0)
      cat("Plotting fits obtained with population predictions\n")
      saemix.plot.fits(x,...)
    },
    "both.fit"={
      x@prefs$level<-c(0:1)
      if(length(x["results"]["ipred"])>0) {
        cat("Plotting the fits overlaying individual and population predictions\n")
        saemix.plot.fits(x,...)
      }
    },
    "residuals.scatter"={
      if(length(x["results"]["wres"])>0) {
        cat("Plotting scatterplots of residuals\n")
        saemix.plot.scatterresiduals(x,...)
      }
    },
    "residuals.distribution"={
      if(length(x["results"]["wres"])>0) {
        cat("Plotting the distribution of residuals\n")
        saemix.plot.distribresiduals(x,...)
      }
    },
    "random.effects"={
      saemix.plot.randeff(x,...)
    },
    "correlations"={
      saemix.plot.correlations(x,...)
    },
    "parameters.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        cat("No covariates in the dataset\n")
        return()
      } else saemix.plot.parcov(x,...)
    },
    "randeff.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        cat("No covariates in the dataset\n")
        return()
      } else saemix.plot.randeffcov(x,...)
    },
    "marginal.distribution"={
      saemix.plot.distpsi(x,...)
    },
    "vpc"={
      if(length(x["results"]["npde"])>0) {
        cat("Plotting VPC\n")
        saemix.plot.vpc(x,...)
      }
    },
    "npde"={
      if(length(x["results"]["npde"])>0) {
        cat("Plotting npde\n")
        saemix.plot.npde(x,...)
      }
    },
    cat("Plot ",ipl," not implemented yet\n")
     )
   }
  }
)

####################################################################################
####			saemixObject class - accesseurs parametres		####
####################################################################################

# Extract individual parameter estimates (psi_i)
setMethod("psi","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
    if(indiv.par=="map") { # mode
      if(length(object@results@map.psi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@map.psi[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.psi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@cond.mean.psi
    }
    return(psi)
  }
)

# Extract individual parameter estimates on non-transformed scale (phi_i)
setMethod("phi","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
    if(indiv.par=="map") { # mode
      if(length(object@results@map.phi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@map.phi[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.phi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@cond.mean.phi
    }
    return(phi)
  }
)

# Extract individual estimates of random effects (eta_i)
setMethod("eta","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
#    cat("Nom objet",namObj,"\n")
    if(indiv.par=="map") { # mode
      if(length(object@results@map.eta)==0) {
        object<-compute.eta.map(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@map.eta[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.eta)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@cond.mean.eta
    }
    return(eta)
  }
)

# Extract coefficients
setMethod("coef","SaemixObject",
  function(object) {
#    if(missing(level)) level<-1
#    if(missing(indiv.par)) indiv.par<-"map"
    pfix<-object@results@fixed.effects[object@results@indx.fix]
    names(pfix)<-object@results@name.fixed[object@results@indx.fix]
#c(object@results@fixed.effects,object@name.res[object@indx.res])
#    names(pfix)<-c(object@results@name.fixed,object@name.res[object@indx.res])
    pop.phi<-object@results@mean.phi
    pop.psi<-transphi(pop.phi,object@model@transform.par)
    ind.psi<-list(map=object@results@map.psi[,-c(1)], cond=object@results@cond.mean.psi)
    ind.phi<-list(map=object@results@map.phi[,-c(1)], cond=object@results@cond.mean.phi)
    eta<-list(map=object@results@map.eta[,-c(1)],cond=object@results@cond.mean.eta)
    colnames(pop.phi)<-colnames(pop.psi)<-colnames(ind.phi$map)<-names(pfix)
    colnames(ind.psi$map)<-colnames(ind.phi$cond)<-colnames(ind.psi$cond)<-names(pfix)
    coef<-list(fixed=pfix,population=list(psi=pop.psi,phi=pop.phi), individual=list(psi=ind.psi,phi=ind.phi,eta=eta))
    return(coef)
  }
)

####################################################################################
