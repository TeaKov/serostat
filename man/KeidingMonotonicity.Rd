\name{Keiding}
\alias{KeidingMonotonicity}
\title{Estimating Force of Infection, According to Keiding}
\usage{
KeidingMonotonicity(p)
}
\arguments{
\item{p}{This is a list containing the following parameters: \bold{pos} describes the amount of successes out of ``tot'' trials, \bold{tot} is the total amount of trials, \bold{grid} gives back the rounded (and unique) age groups, and \bold{bw} that describes the bandwidth. Note that both ``pos'' and ``tot'' should be ordered by age.}
}
\description{
Keiding estimates the force of infection using the PAVA algorithm (see help("pavit")).
}
\examples{
# Using Belgian B19 data
data("VZV_B19_BE_0103")
subset<-(VZV_B19_BE_0103$age>0.5)&(VZV_B19_BE_0103$age<76)&
  (!is.na(VZV_B19_BE_0103$age))&!is.na(VZV_B19_BE_0103$parvores)
VZV_B19_BE_0103<-VZV_B19_BE_0103[subset,]
y<-VZV_B19_BE_0103$parvores[order(VZV_B19_BE_0103$age)]
a<-VZV_B19_BE_0103$age[order(VZV_B19_BE_0103$age)]

grid<-sort(unique(round(a)))
neg<-table(y,round(a))[1,]
pos<-table(y,round(a))[2,]
tot<-neg+pos

params <- list(pos=pos, tot=tot, grid=grid, bw=10)

KeidingMonotonicity(params)

# Using Bulgarian HAV data
data("HAV_BUL_64_BUL_64")
a<-c(rep(HAV_BUL_64$Age,HAV_BUL_64$Pos),rep(HAV_BUL_64$Age,HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-c(rep(rep(1,length(HAV_BUL_64$Age)),HAV_BUL_64$Pos),
  rep(rep(0,length(HAV_BUL_64$Age)),HAV_BUL_64$Tot-HAV_BUL_64$Pos))
y<-y[order(a)]
a<-a[order(a)]

grid<-sort(unique(round(a)))
neg<-HAV_BUL_64$Tot-HAV_BUL_64$Pos
pos<-HAV_BUL_64$Pos
tot<-neg+pos

params <- list(pos=pos, tot=tot, grid=grid, bw=30)

KeidingMonotonicity(params)
}
