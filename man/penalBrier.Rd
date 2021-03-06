% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.r
\name{penalBrier}
\alias{penalBrier}
\title{Evaluation of the performance of risk prediction models with binary status response variable.}
\usage{
penalBrier(Data, coeffP)
}
\arguments{
\item{Data}{a data matrix; in the first column there should be the response variable y.
If you give the training dataset it will calculate the Brier score.}

\item{coeffP}{a named vector of coefficients}
}
\value{
the Brier score (misclassification error)
}
\description{
Evaluation of the performance of risk prediction models with binary status response variable.
}
\details{
Brier score is a measure for classification performance of a binary classifier.
Its values range between [0,1] and the closest is to 0 the better the classifier is.
The area under the curve and the Brier score is used to summarize and compare the performance.
}
\examples{
# use the penalBrier function on a simulated dataset, with given lamda and w.
\dontrun{
set.seed(14)
beta    <- c(3, 2, -1.6, -4)
noise   <- 5
simData <- SimData(N=100,beta=beta, noise=noise, corr=FALSE)

before   <- Sys.time()
stepPenal<- StepPenal(Data=simData, lamda=1.2, w=0.4)
(totRun  <- Sys.time() - before)

(coeff<- stepPenal$coeffP)
 me <- penalBrier(simData,coeff)
}
}
\references{
Brier, G. W. (1950). Verification of forecasts expressed in terms of probability.
Monthly Weather Review 78.
}
