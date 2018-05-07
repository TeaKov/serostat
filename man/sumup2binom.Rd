\name{Convert from Sum-up Data to Binomial Data}
\alias{sumup2binom}
\title{Convert from Sum-up Data to Binomial Data}
\usage{
sumup2binom(pos, ntot)
}
\arguments{
\item{pos}{The amount of successes out of ``tot'' trials.}
\item{ntot}{The total amount of trials.}
}
\description{
This function will take binary data (a vector filled with zeros and ones) and transform it into two numbers: \code{pos} and \code{ntot}, respectfully depicting the amount of positives and the total amount of individuals.

This function will take two numbers: \code{pos} and \code{ntot}, respectfully depicting the amount of positives and the total amount of individuals, and transform it into binary data. The result is a vector of length \code{ntot} with \code{pos} amount of ones and \code{ntot-pos} amount of zeros.
}
\examples{
pos <- 3
ntot <- 10

res <- sumup2binom(pos, ntot)
res
}
