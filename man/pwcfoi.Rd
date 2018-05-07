\name{Piecewise Constant Force of Infection}
\alias{pwcfoi}
\title{Piecewise Constant Force of Infection}
\usage{
pwcfoi(response, x.var, breaks)
}
\arguments{
\item{response}{The response variable.}
\item{x.var}{The predictor variable.}
\item{breaks}{The user chosen categories over which a constant force of infection is assumed.}
}
\description{
Using a piecewise constant force of infection, monotonicity is often achieved by redefining the age categories over which a constant force of infection is assumed.
}
\examples{
# Using Bulgarian HAV data
data("HAV_BUL_64_BUL_64")
a<-c(rep(HAV_BUL_64$Age,HAV_BUL_64$Pos),rep(HAV_BUL_64$Age,HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-c(rep(rep(1,length(HAV_BUL_64$Age)),HAV_BUL_64$Pos),rep(rep(0,length(HAV_BUL_64$Age)),
  HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-y[order(a)]
a<-a[order(a)]

grid<-sort(unique(round(a)))
neg<-HAV_BUL_64$Tot-HAV_BUL_64$Pos
pos<-HAV_BUL_64$Pos
tot<-neg+pos

# monotone curve
breaks1<-c(0,10,20,40,60,100)
pcwfit1<-pwcfoi(response=y,x.var=a,breaks=breaks1)

# non-monotone curve
breaks2<-c(0,5,11,18,60,100)
pcwfit2<-pwcfoi(response=y,x.var=a,breaks=breaks2)
}
