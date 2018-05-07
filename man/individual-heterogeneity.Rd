\name{Individual Heterogeneity - Gamma Failty}
\alias{GF}
\alias{Frailty}
\title{Modeling the Prevalence and the Force of Infection Directly from Antibody Levels}
\usage{
GF(PP, PN, NP, NN, a, alpha1eta, beta1eta,
  gamma1eta, alpha2eta, beta2eta, gamma2eta, thetaeta)
}
\arguments{
\item{PP}{}
\item{PN}{}
\item{NP}{}
\item{NN}{}
\item{a}{}
\item{alpha1eta}{}
\item{beta1eta}{}
\item{gamma1eta}{}
\item{alpha2eta}{}
\item{beta2eta}{}
\item{gamma2eta}{}
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
