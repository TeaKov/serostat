\name{HAV_BUL_64}
\alias{HAV_BUL_64}
\docType{data}
\title{
Hepatitis A Serological Survey Data (Bulgaria, 1964)
}
\description{
Keiding (1991) introduced the Hepatitis A dataset from Bulgaria. The data consists out of a cross-sectional survey conducted in 1964 and contains information about 850 individuals from Bulgaria with age ranging from 1 to 86 years.
}
\usage{data("HAV_BUL_64")}
\format{
  A data frame with 83 observations on the following 3 variables.
  \describe{
    \item{\code{Age}}{A numeric vector containing age.}
    \item{\code{Pos}}{A numeric vector containing the number of seropositives in each age group.}
    \item{\code{Tot}}{A numeric vector containing the total number of individuals in each age group.}
  }
}
\details{
See page 64 of the book.
}
\source{
\url{https://www.jstor.org/stable/2983150}
}
\examples{
data(HAV_BUL_64)
}
\keyword{datasets}
