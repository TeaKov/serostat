\name{Convert from Binomial Data to Sum-up Data}
\alias{binom2sumup}
\title{Convert from Binomial Data to Sum-up Data}
\usage{
binom2sumup(data)
}
\arguments{
\item{data}{The binary data that needs to be summed up.}
}
\description{
This function will take binary data (a vector filled with zeros and ones) and transform it into two numbers: \code{pos} and \code{ntot}, respectfully depicting the amount of positives and the total amount of individuals.
}
\examples{
# Generate binary data
data <- rbinom(10, 1, 0.5)

res <- binom2sumup(data)
res
}
