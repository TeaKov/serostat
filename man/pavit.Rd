\name{PAVIT}
\alias{pavit}
\title{Pool Adjacent Violator Algorithm}
\usage{
pavit(pos,tot)
}
\arguments{
\item{pos}{The amount of successes out of ``tot'' trials. These should be ordered by age.}
\item{tot}{The amount of trials. These should be ordered by age.}
}
\description{
The PAVA algorithm proceeds by recursively checking monotone behavior and by pooling if necessary and finally stops if monotonicity is achieved.
}
\examples{
# Load Belgian B19 data
data("VZV_B19_BE_0103")
subset<-(VZV_B19_BE_0103$age>0.5)&(VZV_B19_BE_0103$age<76)&
  (!is.na(VZV_B19_BE_0103$age))&!is.na(VZV_B19_BE_0103$parvores)
VZV_B19_BE_0103<-VZV_B19_BE_0103[subset,]
y<-VZV_B19_BE_0103$parvores[order(VZV_B19_BE_0103$age)]
a<-VZV_B19_BE_0103$age[order(VZV_B19_BE_0103$age)]

neg<-table(y,round(a))[1,]
pos<-table(y,round(a))[2,]
tot<-neg+pos

pavit(pos=pos, tot=tot)
}
