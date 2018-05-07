\name{Semiparametric Regression - P-Splines}
\alias{pspline.fit}
\alias{pspline.fitter}
\alias{pspline.predictor}
\alias{pspline.checker}
\title{Semiparametric Approaches to Model the Prevalence and Force of Infection}
\usage{
pspline.fit(response, x.var, ps.intervals, wts, degree, order,
	link, family, m.binomial, r.gamma, lambda, x.predicted, ridge.adj)
}
\arguments{
\item{response}{Response variable}
\item{x.var}{Explanatory variable on abcissae}
\item{ps.intervals}{Number of intervals for B-splines. Default=8.}
\item{wts}{Vector of weights; default is vector of ones.}
\item{degree}{Degree of B-splines. Default=3.}
\item{order}{Order of difference penalty. Default=3.}
\item{link}{Link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical).}
\item{family}{What kind of distribution (family=gaussian, binomial, poisson, gamma)}
\item{m.binomial}{Vector of binomial trials. Default is 1 vector.}
\item{r.gamma}{Vector of gamma shape parameters. Default is 1 vector.}
\item{lambda}{Smoothness regularizing parameter ( >= 0). Default=0.}
\item{x.predicted}{A list of x variables for prediction and twice stderr limits.}
\item{ridge.adj}{Default=0.0001}
}
\description{
Univariate smoother using P-splines, by Paul Eilers & Brian Marx (1995). For the smooth splines, thin plate splines, and adaptive smoothing splines approaches, go to help("SemiParametricRegression").
}
\examples{
# Load Belgian B19 data
data("VZV_B19_BE_0103")
subset<-(VZV_B19_BE_0103$age>0.5)&(VZV_B19_BE_0103$age<76)&(!is.na(VZV_B19_BE_0103$age))&
  !is.na(VZV_B19_BE_0103$parvores)
VZV_B19_BE_0103<-VZV_B19_BE_0103[subset,]
y<-VZV_B19_BE_0103$parvores[order(VZV_B19_BE_0103$age)]
a<-VZV_B19_BE_0103$age[order(VZV_B19_BE_0103$age)]
s<-VZV_B19_BE_0103$sex[order(VZV_B19_BE_0103$age)]

grid<-sort(unique(round(a)))
neg<-table(y,round(a))[1,]
pos<-table(y,round(a))[2,]
tot<-neg+pos

degree <- 3
order <- 2

# P-splines, with logit link-function
lopt <- 80
pspline.fit(response=y, x.var=a, ps.intervals=20, degree=degree, order=order,
  link="logit", family="binomial", lambda=lopt, x.predicted=a)

# P-splines, with cloglog link-function
lopt <- 200
pspline.fit(response=y, x.var=a, ps.intervals=20, degree=degree, order=order,
  link="cloglog", family="binomial", lambda=lopt, x.predicted=a)
}
