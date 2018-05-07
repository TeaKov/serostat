\name{Piecewise Constant Rate}
\alias{pwcrate}
\title{Piecewise Constant Rate}
\usage{
pwcrate(y.var, x.var, n.var, breaks, startpar)
}
\arguments{
\item{y.var}{The response variable (binary, 1=past infection, 0=otherwise).}
\item{x.var}{Numeric vector of age categories.}
\item{n.var}{Number of individuals.}
\item{breaks}{Age categories to be considered.}
\item{startpar}{Starting values for the \dQuote{nlm} method used in this function. This is a vector of length(breaks)-1.}
}
\description{
Using a piecewise constant force of infection, monotonicity is often achieved by redefining the age categories over which a constant force of infection is assumed.
}
\examples{
# Using Belgian B19 data
data("VZV_B19_BE_0103")
VZV_B19_BE_0103 <- VZV_B19_BE_0103[!is.na(VZV_B19_BE_0103$parvores)&!is.na(VZV_B19_BE_0103$age)&
  VZV_B19_BE_0103$age<70&VZV_B19_BE_0103$age>=1,]
VZV_B19_BE_0103 <- VZV_B19_BE_0103[order(VZV_B19_BE_0103$age),]

y <- VZV_B19_BE_0103$parvores
a <- VZV_B19_BE_0103$age

# Age category
breakpoints <- c(0.5,2,6,12,19,31,100)

# Monotonized piecewise constant FOI
pcw.fit <- pwcrate(y.var=y, x.var=a, breaks=breakpoints)
}
