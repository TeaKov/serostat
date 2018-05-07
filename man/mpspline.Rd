\name{Penalized B-splines with Monotonicity Constraint}
\alias{mpspline.fit}
\alias{Threshold}
\title{Penalized B-splines with Monotonicity Constraint}
\usage{
mpspline.fit(response, x.var, ps.intervals, degree, order, link, family, alpha, kappa)
}
\arguments{
\item{response}{Response variable}
\item{x.var}{Explanatory variable on abcissae}
\item{ps.intervals}{Number of intervals for B-splines. Default=20.}
\item{degree}{Degree of B-splines. Default=3.}
\item{order}{Order of difference penalty. Default=2.}
\item{link}{Link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical). Default=logit}
\item{family}{What kind of distribution (family=gaussian, binomial, poisson, gamma). Default=binomial.}
\item{alpha}{Smoothness regularizing parameter (>= 0). Default=2.}
\item{kappa}{Imposes to what extent non-monotone behaviour is penalized. Typically chosen to be large. Default=1e8.}
}
\description{
This is an adjusted version of the original ``pspline.fit''-function, taking monotonicity constraints into consideration (Eilers et al., 1996). For additional documentation regarding the smoothing splines, thin plate splines, and adaptive smoothing splines, go to \code{\link{SemiParametricRegression}}).
}
\examples{
# Load Bulgarian HAV data
data(HAV_BUL_64)
a<-c(rep(HAV_BUL_64$Age,HAV_BUL_64$Pos),
  rep(HAV_BUL_64$Age,HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-c(rep(rep(1,length(HAV_BUL_64$Age)),HAV_BUL_64$Pos),
  rep(rep(0,length(HAV_BUL_64$Age)),HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-y[order(a)]
a<-a[order(a)]
grid<-sort(HAV_BUL_64$Age)
neg<-HAV_BUL_64$Tot-HAV_BUL_64$Pos
pos<-HAV_BUL_64$Pos
tot<-neg+pos

# P-splines, with logit link-function
fit0<-mpspline.fit(response=y, x.var=a, ps.intervals=20,
  degree=3, order=2, link="logit", family="binomial", alpha=54, kappa=0)
fit1<-mpspline.fit(response=y, x.var=a, ps.intervals=20,
  degree=3, order=2, link="logit", family="binomial", alpha=66, kappa=1e8)

\donttest{dev.new(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(2, 0.5, 0),mar=c(4.1,4.1,4.1,3))
plot(grid,pos/tot,cex=0.05*tot,pch=1,xlab="age",
  ylab="seroprevalence",xlim=c(0,72),ylim=c(0,1),lwd=2)
lines(fit0$x,fit0$y,lwd=3,lty=1)
lines(fit1$x,fit1$y,lwd=3,lty=3)
lines(foi.num(fit1$x,fit1$y)$grid,foi.num(fit1$x,fit1$y)$foi,lwd=3,lty=1)
lines(foi.num(fit0$x,fit0$y)$grid,foi.num(fit0$x,fit0$y)$foi,lwd=3,lty=3)}

# Plotting logarithm of antibody activity levels in U/ml as a function of
# the individual's age, with threshold values represented by solid horizontal lines.
# Upped dashed horizontal line: estimated mean for the infected population.
# Lower dashed horizontal line: estimated mean for the susceptible population.
# Solid smooth curve: monotone least-squares fit using p-splines
#
# Using Belgian B19 data
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age<40.5)&(!is.na(VZV_B19_BE_0103$age))&
  (!is.na(VZV_B19_BE_0103$VZVmUIml))
data <- VZV_B19_BE_0103[subset,]

# Data to use when taking the continuous levels
z<-log(data$VZVmUIml[order(data$age)]+1)
a<-data$age[order(data$age)]

\donttest{dev.new(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(3, 0.5, 0))

plot(a,z,main="",pch=16,xlab="age",ylab="log(antibody levels+1)",
  xlim=c(0,40.5),cex=0.5,col="dark grey")
abline(h=log(51))
abline(h=log(101))
abline(h=2.316,lty=5,lwd=2)  #mean of lower component of mixture, fitted in code below
abline(h=6.338,lty=5,lwd=2)}

# Selecting the best smoothing parameter with BIC
alphagrid <- seq(0.01,1,by=0.01)
res <- matrix(ncol=2,nrow=length(alphagrid))
for (i in (1:length(alphagrid))) {
  fitC <- mpspline.fit(response=z, x.var=a, ps.intervals=20, degree=3,
    link="identity", family="gaussian", order=2, alpha=alphagrid[i], kappa=1e8)
  res[i,1] <- alphagrid[i]
  res[i,2] <- fitC$bic
}
alphafin <- res[res[,2] == min(res[,2]),1]
fitC <- mpspline.fit(response=z, x.var=a, ps.intervals=20, degree=3,
  link="identity", family="gaussian", order=2, alpha=alphafin, kappa=1e8)
\donttest{lines(fitC$x, fitC$y, lwd=3, lty=1)}

# Proportion positive, as a function of the corresponding half-year
# age categories, overlaid with the monotone p-spline fit and the FOI:
# with BIC optimal smoothing parameter
#
# Using Belgian B19 data
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age<40.5)&(!is.na(VZV_B19_BE_0103$age))&
  (!is.na(VZV_B19_BE_0103$VZVmUIml)&(!is.na(VZV_B19_BE_0103$VZVres)))
data <- VZV_B19_BE_0103[subset,]

# Data to use when taking the binary indicators (different because some are inconclusive)
y <- data$VZVres[order(data$age)]
a <- data$age[order(data$age)]

alphagrid <- seq(1,100,by=1)
res <- matrix(ncol=2, nrow=length(alphagrid))
for (i in (1:length(alphagrid))) {
  fit1 <- mpspline.fit(response=y, x.var=a, ps.intervals=20, degree=3, order=2,
    link="logit", family="binomial", alpha=alphagrid[i], kappa=1e8)
  res[i,1] <- alphagrid[i]
  res[i,2] <- fit1$bic
}

alphafin <- res[res[,2] == min(res[,2]),1]
fit1 <- mpspline.fit(response=y, x.var=a, ps.intervals=20, degree=3, order=2,
  link="logit", family="binomial", alpha=alphafin, kappa=1e8)
grid<-sort(unique(round(a)))
neg<-table(y,round(a))[1,]
pos<-table(y,round(a))[2,]
tot<-neg+pos

\donttest{dev.new(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(2.5, 0.5, 0),mar=c(5.1,3.5,4.1,4))

plot(grid,pos/tot,cex=0.02*tot,xlab="age",ylab="seroprevalence",
  xlim=c(0,max(grid)),ylim=c(0,1),lwd=2)
lines(fit1$x,fit1$y,lwd=2,lty=1)
lines(foi.num(fit1$x,fit1$y)$grid,foi.num(fit1$x,fit1$y)$foi,lwd=2,lty=1)
axis(side=4,at=c(0.0,0.1,0.2,0.3,0.4))
mtext(side=4,las=3,"force-of-infection",line=1.5)}
}
\seealso{\cite{Eilers, P. H. C.; Marx, B. D. (1996). Flexible smoothing with B-splines and penalties, 11(2), 89-121.}}
