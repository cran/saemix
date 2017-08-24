\docType{methods}
\name{plot-method}
\alias{plot,ANY-method}
\alias{plot.saemix}
\alias{plot,SaemixData}
\alias{plot,SaemixData-method}
\alias{plot,SaemixSimData}
\alias{plot,SaemixSimData-method}
\title{Methods for Function plot}
\description{
Methods for function \code{plot}
}
\usage{
\S4method{plot}{SaemixData}(x, y, ...)
}
\arguments{
\item{x}{an object of class SaemixData, SaemixSimData or SaemixObject}

\item{y}{empty}

\item{\dots}{optional arguments passed to the plots}
}
\section{Methods}{

\describe{
\item{list("signature(x = \"ANY\")")}{ default plot function ?}
\item{list("signature(x = \"SaemixData\")")}{ Plots the data. Defaults to a
spaghetti plot of response versus predictor, with lines joining the data for
one individual.}
\item{list("signature(x = \"SaemixObject\")")}{ This method gives access to
a number of plots that can be performed on a SaemixObject}
\item{list("signature(x = \"SaemixSimData\")")}{ Plots simulated datasets} }

}
\keyword{methods}

