% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.r
\name{tuneParamCL2}
\alias{tuneParamCL2}
\title{Tune parameters w and lamda using the CL2 penalty}
\usage{
tuneParamCL2(Data, nfolds = nfolds, grid, algorithm = c("QN"))
}
\arguments{
\item{Data}{a data frame, as a first column should have the response variable y and the other columns the predictors}

\item{nfolds}{the number of folds used for cross-validation. OBS! nfolds>=2}

\item{grid}{a grid (data frame) with values of lamda and w that will be used for tuning
to tune the model. It is created by expand.grid see example below}

\item{algorithm}{choose between BFGS ("QN") and hjk (Hooke-Jeeves optimization free) to be used for optmization}
}
\value{
A matrix with the following: the average (over folds) cross-validated AUC, the totalVariables selected on the training set,
and the standard deviation of the AUC over the nfolds
}
\description{
Does k-fold cross-validation with the function optimPenalLikL2 and returns
the values of lamda and w that maximize the area under the ROC.
}
\details{
It supports the BFGS optimization method ('QN') from the optim {stats} function, the Hooke-Jeeves derivative-free
minimization algorithm ('hjk').
The value of lamda and w that yield the maximum AUC on the
cross-validating data set is selected. If more that one value of lamda nad w yield the same AUC,
then the biggest values of lamda and w are choosen.
}
