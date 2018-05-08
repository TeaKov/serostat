\name{Individual Heterogeneity - Gamma Failty}
\alias{GF}
\alias{Frailty}
\alias{GammaFrailty}
\alias{IndividualHeterogeneity}
\title{Modeling the Prevalence and the Force of Infection Directly from Antibody Levels}
\usage{
GF(PP, PN, NP, NN, a, alpha1eta, beta1eta,
  gamma1eta, alpha2eta, beta2eta, gamma2eta, thetaeta)
}
\arguments{
\item{PP}{A numeric vector containing the number of seropositives for both diseases.}
\item{PN}{A numeric vector containing the number of individuals that are seropositive for the first disease and seronegative for the second disease.}
\item{NP}{A numeric vector containing the number of individuals that are seropositive for the second disease and seronegative for the first disease.}
\item{NN}{A numeric vector containing the number of seronegatives for both diseases.}
\item{a}{A numeric vector containing age.}
\item{alpha1eta}{Nonnegative scaling parameter for first disease.}
\item{beta1eta}{Nonnegative scaling parameter for first disease.}
\item{gamma1eta}{Long-term residual value of the force of infection for the first disease (nonnegative).}
\item{alpha2eta}{Nonnegative scaling parameter for second disease.}
\item{beta2eta}{Nonnegative scaling parameter for second disease.}
\item{gamma2eta}{Long-term residual value of the force of infection for the second disease (nonnegative).}
\item{thetaeta}{The shape parameter of the gamma frailty distribution.}
}
\description{
Individual heterogeneity has been shown to be a key factor in the estimation of the basic reproduction number. Individuals are dissimilar in the way they acquire infections. Some individuals are more susceptible than others and will experience infection earlier.
}
\examples{
# Load Belgian B19 data
library(stats4)
data("VZV_B19_BE_0103")
VZV_B19_BE_0103 <- VZV_B19_BE_0103[!is.na(VZV_B19_BE_0103$VZVres)&
  !is.na(VZV_B19_BE_0103$parvores)&!is.na(VZV_B19_BE_0103$age)&
  VZV_B19_BE_0103$age<70&VZV_B19_BE_0103$age>=1,]
VZV_B19_BE_0103 <- VZV_B19_BE_0103[order(VZV_B19_BE_0103$age),]

y1 <- VZV_B19_BE_0103$VZVres
y2 <- VZV_B19_BE_0103$parvores
age <- VZV_B19_BE_0103$age

a <- unique(age)
covariate <- seq(min(age),max(age),1)

# Counts per age-value
PP <- as.vector(hist(age[y1=="1"&y2=="1"], plot=FALSE, breaks=c(0,a))$counts)
PN <- as.vector(hist(age[y1=="1"&y2=="0"], plot=FALSE, breaks=c(0,a))$counts)
NP <- as.vector(hist(age[y1=="0"&y2=="1"], plot=FALSE, breaks=c(0,a))$counts)
NN <- as.vector(hist(age[y1=="0"&y2=="0"], plot=FALSE, breaks=c(0,a))$counts)

\donttest{fit <- mle(GF, start=list(alpha1eta=-1, beta1eta=-1,
  gamma1eta=-1, alpha2eta=-1, beta2eta=-1, gamma2eta=-1, thetaeta=-1),
  fixed=list(PP=PP, PN=PN, NP=NP, NN=NN, a=a))
summary(fit)}
}
\details{
See page 186 of the book.
}
\seealso{\code{\link{vgam}}}
