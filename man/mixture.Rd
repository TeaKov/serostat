\name{Mixtures}
\alias{Mixtures}
\alias{mix}
\title{Modeling the Prevalence and the Force of Infection Directly from Antibody Levels}
\usage{
mix(mixdat, mixpar, dist = "norm", \ldots)
}
\arguments{
\item{mixdat}{A data frame containing grouped data, whose first column should be right boundaries of grouping intervals where the first and last intervals are open-ended; whose second column should consist of the frequencies indicating numbers of observations falling into each interval. If conditional data are available, this data frame should have k + 2 columns, where k is the number of components, whose element in row j and column i + 2 is the number of observations from the jth interval belonging to the ith component.}
\item{mixpar}{A data frame containing starting values for parameters of component distributions, which are, in order, the proportions, means, and standard deviations.}
\item{dist}{The distribution of components, it can be one of \dQuote{norm}, \dQuote{lnorm}, \dQuote{gamma}, \dQuote{weibull}, \dQuote{binom}, \dQuote{nbinom}, and \dQuote{pois}}
\item{\ldots}{See \code{\link{mix}}.}
}
\description{
True prevalence and force of infection can be derived directly by using a two-component mixture model for the antibody activity levels. The group of infected individuals and the group of susceptible individuals form the two components with the true prevalence as age-dependent mixing probability.
}
\examples{
library(mixdist)
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age<40.5)&(!is.na(VZV_B19_BE_0103$age))&
  (!is.na(VZV_B19_BE_0103$VZVmUIml))
VZV_B19_BE_0103<-VZV_B19_BE_0103[subset,]

z1<-log(VZV_B19_BE_0103$VZVmUIml[order(VZV_B19_BE_0103$age)]+1)

zmixdat <- mixgroup(z1,breaks=40)
zstartpar <- mixparam(pi=c(0.2,0.8), mu=c(2,6), sigma=c(0.5,1))
mixfit <- mix(zmixdat,zstartpar,dist="norm")
\donttest{summary(mixfit)}
}
\details{
See page 171 of the book.
}
\seealso{\code{\link{mpspline.fit}}, \code{\link{Threshold}}, \code{\link{mix}}}
