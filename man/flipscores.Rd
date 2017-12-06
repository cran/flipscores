\name{flipscores}
\alias{flipscores}
\title{Robust testing in GLMs}
\description{
  Provides two robust tests for testing in GLMs, by sign-flipping score contributions. The tests are often robust against overdispersion, heteroscedasticity and, in some cases, ignored nuisance variables.
}
\usage{
  flipscores(model0, model1, X1, alternative = "two.sided",  w=1E5, scoretype="basic")
}
\arguments{
  \item{model0}{A glm or lm object representing the fit under H0.}
  \item{model1}{A glm or lm object representing the fit under H1.
   It contains one covariate in addition to those in model0.
   If that covariate is a factor, it should have at most 2 levels.
   Either \code{model1} or \code{X1} should be specified.}
  \item{X1}{Vector with the values of the covariate of interest.
  If X1 is a factor, it should have at most 2 levels.
  Either \code{model1} or \code{X1} should be specified.}
  \item{alternative}{Should be \code{"greater"}, \code{"less"} or \code{"two.sided"}}
  \item{w}{The number of times that the scores are randomly sign-flipped.}
  \item{scoretype}{The type of score that is computed, either \code{"basic"} or \code{"effective"}.
  Using \code{"effective"} takes into account nuisance estimation.}
}

\value{A p-value.}

\examples{

set.seed(8153)
n <-50
beta <- 2  #coeffient of interest (H0: beta=0)
gamma <- 1   #nuisance coefficient

dataset <- data.frame(x=NULL,z=NULL,y=NULL)
dataset[n,] <- NA
dataset[,1:2] <- data.frame(matrix(rnorm(n*2),n,2))  #generate covariates
colnames(dataset) <- c("x","z")
dataset$y <- rnbinom(n, mu = exp(dataset[,1]*beta + dataset[,2]*gamma), size=1 )

#Y has a negative binomial distribution but we assume a Poisson model:
modelz <- glm(y~z, family=poisson,data=dataset,x=TRUE)


## Basic test:
pv1 <- flipscores(model0=modelz, X1=dataset[,1], alternative = "two.sided", scoretype="basic")
pv1  #p-value

## Test that takes into account nuisance estimation:
pv2 <- flipscores(model0=modelz, X1=dataset[,1], alternative = "two.sided", scoretype="eff")
pv2  #p-value

}

