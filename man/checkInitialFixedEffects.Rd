% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SaemixModel.R
\name{checkInitialFixedEffects}
\alias{checkInitialFixedEffects}
\title{Check initial fixed effects for an SaemixModel object applied to an SaemixData object}
\usage{
checkInitialFixedEffects(model, data, psi = c(), id = c(), ...)
}
\arguments{
\item{model}{an SaemixModel object}

\item{data}{an SaemixData object (the predictors will then be extracted from the object using the name.predictors slot of the object)}

\item{psi}{a vector or a dataframe giving the parameters for which predictions are to be computed (defaults to empty).
The number of columns in psi (or the number of elements of psi, if psi is given as a vector) should match the number of
parameters in the model, otherwise an error message will be shown and the function will return empty.
If psi is NA, the predictions are computed for the population parameters in the model (first line of the psi0 slot).
Covariates are not taken into account in the prediction.
If psi is a dataframe, each line will be used for a separate 'subject' in the predictors dataframe, as
indicated by the id argument; if id is not given, only the first line of psi will be used.}

\item{id}{the vector of subjects for which individual plots will be obtained. If empty, the first 12 subjects in the dataset will be used (subject id's are taken from the name.group slot in the data object). If id is given, individual plots will be shown for the matching subjects in the dataset (eg if id=c(1:6), the first 6 subjects in the dataframe will be used for the plots, retrieving their ID from the data object)}

\item{\dots}{unused argument, for consistency with the generic}
}
\value{
the predictions corresponding to the values for each observation in the predictors of either the model f or log-likelihood.
For Gaussian data models, the function also plots the data overlayed with the model predictions for each subject in id
(where id is the index in the N subjects).
}
\description{
Check initial fixed effects for an SaemixModel object applied to an SaemixData object
}
\details{
The function uses the model slot of the SaemixModel object to obtain predictions, using the predictors object. The
user is responsible for giving all the predictors needed by the model function.
if psi is not given, the predictions will be computed for the population parameters (first line of the psi0 slot) of the object.

The predictions correspond to the structure of the model. For models defined as a structural model,
individual plots for the subjects in id overlaying the predictions for the parameters psi and the individual data
are shown, and the predictions correspond to f(t_ij, psi).
For models defined in terms of their likelihood, the predictions returned correspond to the log-likelihood.
No individual graphs are currently available for discrete data models.

Warning: this function is currently under development and the output may change in future versions of the package
}
\examples{
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
	  dose<-xidep[,1]
	  tim<-xidep[,2]  
	  ka<-psi[id,1]
	  V<-psi[id,2]
	  CL<-psi[id,3]
	  k<-CL/V
	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	  return(ypred)
}

saemix.model<-saemixModel(model=model1cpt,modeltype="structural",
  description="One-compartment model with first-order absorption", 
  psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
  dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
  covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6))
checkInitialFixedEffects(saemix.model, saemix.data, id=c(1:6), psi=c(0.5, 30, 2)) # better fit

}
