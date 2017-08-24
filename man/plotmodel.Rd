\name{plotmodel}
\alias{plotmodel,SaemixModel-method}
\alias{plotmodel,SaemixModel}
\alias{plotmodel}
\title{Plot output from a SaemixModel}
\usage{
plotmodel(x, y, range = c(0, 1), psi, predictors, ...)
}
\arguments{
\item{x}{an SaemixModel object, for example a model returned by the \code{\link{saemixModel}} function or the model component of a SaemixObject object}

\item{y}{empty}

\item{range}{range for the predictor shown on the X-axis of the plot}

\item{psi}{value of the model parameters}

\item{predictors}{value of the other predictors in the model}

\item{\dots}{optional arguments passed to the plots}
}
\description{
Plots the predictions from a model.
}
\value{
Plot of the model predictions versus the main predictor in the model, given values for the model parameters and additional predictors.
}
\keyword{plots,methods}

