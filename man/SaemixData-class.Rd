\name{SaemixData-class}
\Rdversion{1.1}
\docType{class}
\alias{SaemixData-class}
\alias{SaemixData}
\alias{[<-,SaemixData-method}
\alias{[,SaemixData-method}
\alias{initialize,SaemixData-method}
\alias{plot,SaemixData-method}
\alias{print,SaemixData-method}
\alias{read.saemixData,SaemixData-method}
\alias{showall,SaemixData-method}
\alias{show,SaemixData-method}
\alias{plot,SaemixData}
\alias{print,SaemixData}
\alias{showall,SaemixData}
\alias{show,SaemixData}
\alias{summary,SaemixData}
\alias{summary,SaemixData-method}

\alias{SaemixRepData-class}
\alias{SaemixRepData}
\alias{[<-,SaemixRepData-method}
\alias{[,SaemixRepData-method}
\alias{initialize,SaemixRepData-method}
\alias{show,SaemixRepData-method}

\alias{SaemixSimData-class}
\alias{SaemixSimData}
\alias{[<-,SaemixSimData-method}
\alias{[,SaemixSimData-method}
\alias{initialize,SaemixSimData-method}
\alias{plot,SaemixSimData-method}
\alias{show,SaemixSimData-method}

\title{Class "SaemixData"}
\description{An object of the SaemixData class, representing a longitudinal data structure, used by the SAEM algorithm.
}
\section{Objects from the Class}{An object of the SaemixData class can be created by using the function \code{\link{saemixData}}
}
\section{Slots}{
  \describe{
    \item{\code{name.data}:}{Object of class \code{"character"}: name of the dataset}
    \item{\code{header}:}{Object of class \code{"logical"}: whether the dataset/file contains a header. Defaults to TRUE }
    \item{\code{sep}:}{Object of class \code{"character"}: the field separator character}
    \item{\code{na}:}{Object of class \code{"character"}: a character vector of the strings which are to be interpreted as NA values}
    \item{\code{name.group}:}{Object of class \code{"character"}: name of the column containing the subject id}
    \item{\code{name.predictors}:}{Object of class \code{"character"}: name of the column(s) containing the predictors}
    \item{\code{name.response}:}{Object of class \code{"character"}: name of the column containing the response variable y modelled by predictor(s) x}
    \item{\code{name.covariates}:}{Object of class \code{"character"}: name of the column(s) containing the covariates, if present (otherwise empty)}
    \item{\code{name.X}:}{Object of class \code{"character"}: name of the column containing the regression variable to be used on the X axis in the plots}
    \item{\code{units}:}{Object of class \code{"list"}: list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates}
    \item{\code{tab}:}{Object of class \code{"data.frame"}: dataframe containing the data}
    \item{\code{N}:}{Object of class \code{"numeric"}: number of subjects}
    \item{\code{index}:}{Object of class \code{"numeric"}: subject index}
    \item{\code{id}:}{Object of class \code{"character"}: subject id data}
    \item{\code{xind}:}{Object of class \code{"data.frame"}: predictors data (eg dose, time)}
    \item{\code{y}:}{Object of class \code{"numeric"}: response data}
    \item{\code{yorig}:}{Object of class \code{"numeric"}: response data, on the original scale}
    \item{\code{cov}:}{Object of class \code{"data.frame"}: covariates, if present in the model (otherwise empty)}
    \item{\code{ntot.obs}:}{Object of class \code{"numeric"}: total number of observations}
    \item{\code{nind.obs}:}{Object of class \code{"numeric"}: vector containing the number of observations for each subject}
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "SaemixData")}: replace elements of object}
    \item{[}{\code{signature(x = "SaemixData")}: access elements of object}
    \item{initialize}{\code{signature(.Object = "SaemixData")}: internal function to initialise object, not to be used}
    \item{plot}{\code{signature(x = "SaemixData")}: plot the data}
    \item{print}{\code{signature(x = "SaemixData")}: details about the object}
    \item{read.saemixData}{\code{signature(object = "SaemixData")}: internal function, not to be used }
    \item{showall}{\code{signature(object = "SaemixData")}: all the elements in the object}
    \item{show}{\code{signature(object = "SaemixData")}: details about the object}
    \item{summary}{\code{signature(object = "SaemixData")}: summary of the data}
	 }
}
\references{ Kuhn, E., and Lavielle, M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Monolix32_UsersGuide.pdf (\url{http://software.monolix.org/sdoms/software/})
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{saemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}}
\examples{
showClass("SaemixData")

# Specifying column names
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# Specifying column numbers
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
  name.group=1,name.predictors=c(2,3),name.response=c(4), name.covariates=5:6, 
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# No column names specified, using automatic recognition of column names
data(PD1.saemix)
saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE, 
  name.covariates=c("gender"),units=list(x="mg",y="-",covariates=c("-")))
}
\keyword{classes}
