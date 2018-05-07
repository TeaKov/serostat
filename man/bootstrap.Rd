\name{Bootstrapping}
\alias{bootstrapping}
\alias{boot}
\title{Bootstrapping}
\usage{
boot(data, statistic, R, \ldots)
}
\arguments{
\item{data}{The data as a vector, matrix or data frame. If it is a matrix or data frame then each row is considered as one multivariate observation.}
\item{statistic}{A function which when applied to data returns a vector containing the statistic(s) of interest. When \code{sim = "parametric"}, the first argument to \code{statistic} must be the data. For each replicate a simulated dataset returned by \code{ran.gen} (see \code{\link{boot}}) will be passed. In all other cases \code{statistic} must take at least two arguments. The first argument passed will always be the original data. The second will be a vector of indices, frequencies or weights which define the bootstrap sample. Further, if predictions are required, then a third argument is required which would be a vector of the random indices used to generate the bootstrap predictions. Any further arguments can be passed to \code{statistic} through the \code{\ldots} argument.}
\item{R}{The number of bootstrap replicates. Usually this will be a single positive integer. For importance resampling, some resamples may use one set of weights and others use a different set of weights. In this case R would be a vector of integers where each component gives the number of resamples from each of the rows of weights.}
\item{\ldots}{See \code{\link{boot}}.}
}
\description{
Bootstrap methods are simulation based resampling techniques for assessing distributional properties of an estimator, such as bias, variability, quantiles, and percentiles. They are particularly useful when standard maximum likelihood inference is complex or unavailable, but they can also be applied to verify standard approximations.

In the nonparametric bootstrap, we generate bootstrap data from the empirical distribution function (EDF). This estimate, which puts mass \code{1/n} on each of the observed \eqn{latex}{y_i}, is a nonparametric estimator for the unknown distribution function F. It is known to be a consistent estimator for F (under some mild regularity conditions). Generating data from the EDF is nothing else than sampling from the original data with replacement, so it is resampling the original sample (a bit as recycling).

In contrast to the nonparametric bootstrap, one has to assume a particular parametric distribution (with estimated parameters) to generate the bootstrap data. The choice of an appropriate distribution is of course crucial, and misspecification of the generating distribution might lead to substantial bias.
}
\examples{
# Load Belgian VZV data
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age>0.5)&(VZV_B19_BE_0103$age<76)&
  (!is.na(VZV_B19_BE_0103$age))&(!is.na(VZV_B19_BE_0103$VZVres))&
  (VZV_B19_BE_0103$VZVmUIml>0.)
VZV_B19_BE_0103 <- VZV_B19_BE_0103[subset,]
d <- VZV_B19_BE_0103$VZVres[order(VZV_B19_BE_0103$age)]
z <- log(VZV_B19_BE_0103$VZVmUIml[order(VZV_B19_BE_0103$age)])
a <- VZV_B19_BE_0103$age[order(VZV_B19_BE_0103$age)]

library(boot)

########
# Nonparametric bootstrapping
########
# For the mean antibody level
meanstat <- function(data,indices) {
  data <- data[indices]
  mean(data)
}

mub=boot(z,meanstat,R=999)
boot.ci(mub,type = c("norm","basic","perc"))

# For the prevalence
pib=boot(d,meanstat,R=999)
boot.ci(pib,type = c("norm","basic","perc"))

# For the effect of age
slope <- function(data,indices) {
  data <- data[indices,]
  d <- data[,2]
  age <- data[,1]
  fit <- glm(d~age,family=binomial)
  fit$coef[2]
}

data <- data.frame(cbind(a,d))
slopeb <- boot(data,slope,R=999)
boot.ci(slopeb,type = c("norm","basic","perc"))

########
# Parametric bootstrapping
########
# Function to generate normal data; mle will contain
# the mean and standard deviation of the original data
z.rg1 <- function(data,mle) {
  out <- data
  out <- rnorm(length(out),mle[[1]],mle[[2]])
  out
}

mub <- boot(z, meanstat, sim="parametric", ran.gen=z.rg1,
  mle=list(mn=mean(z), sd=sqrt(var(z))), R=999)
boot.ci(mub,type = c("norm","basic","perc"))

# Function to generate uniform data, obviously a bad choice!
z.rg2 <- function(data,mle) {
  out <- data
  out <- runif(length(out),0,mle)
  out
}

mub<- boot(z, meanstat, sim="parametric", ran.gen=z.rg2,
  mle=mean(z), R=999)
boot.ci(mub,type = c("norm","basic","perc"))
}
\seealso{\code{\link{boot}}, \code{\link{boot.ci}}}
