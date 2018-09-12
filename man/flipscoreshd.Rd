\name{flipscoreshd}
\alias{flipscoreshd}
\title{Robust testing in GLMs}
\description{
  Similar to the flipscores() function, but allows a multi-dimensional (even high-dimensional)
  parameter of interest.
}
\usage{
  flipscoreshd(model0, X1,  w=200, scoretype="basic", V="identity")
}
\arguments{
  \item{model0}{A glm or lm object representing the fit under H0.}
  \item{X1}{n by d dataframe or matrix with the d tested covariates (which have no effect under H0).}
  \item{w}{The number of times that the scores are randomly sign-flipped.}
  \item{scoretype}{The type of score that is used, either \code{"basic"} or \code{"effective"}.
  Using \code{"effective"} takes into account nuisance estimation.}
  \item{V}{User-defined d-by-d matrix, where d is the dimension of the parameter of
  interest, i.e. the nr of columns of X1.
  One can also specify \code{V="identity"} to use the identity matrix,
  or specify \code{V="invinfo"} to use the inverse of the effective Fisher info.}
}

\value{A p-value.}

\examples{

set.seed(8193)
n <-50

dataset= data.frame(matrix(rnorm(n*5),n,5))
  
for(k in 1:n){
  dataset[k,1:5] <- dataset[k,1:5] + rnorm(1)   #make the covariates correlated
}

#Y has a negative binomial distribution but we assume a Poisson model:
#Y depends only on the 4-th covariate:
dataset$y <- rnbinom(n, mu = exp(dataset[,4]*0.5), size=1 )

#We will test the hypothesis H0 that X1,X2,X3 all have coefficient 0.
#So H0 is true. X4,X5 are nuisance covariates.

formula0 <- as.formula("y ~ X4+X5")           #X4,X5: nuisance 
X1 <- dataset[,1:3]                           #X1,X2,X3: tested

model0 = glm(formula0, family=poisson,na.action=na.omit,data=dataset,x=TRUE)

## Basic test (often conservative):
pv1 <- flipscoreshd(model0=model0, X1=X1, w=200, scoretype="basic",V="invinfo")
pv1  #p-value

## Test that takes into account nuisance estimation:
pv1 <- flipscoreshd(model0=model0, X1=X1, w=200, scoretype="effective",V="invinfo")
pv1  #p-value

}

