\name{HAV_BE_02}
\alias{HAV_BE_02}
\docType{data}
\title{
Hepatitis A Serological Survey Data (Flanders, Belgium, 2002)
}
\description{
After the serological study conducted between 1993 and early 1994 (the HAV_BE_9394 dataset), a second serological survey on Hepatitis A seroprevalence has been conducted in 2002. This serological sample is a subset of the serological dataset on Varicella-Zoster Virus and B19 serology (ESEN2 study) in Belgium (VZV_B19_BE_0103 dataset) where only individuals living in Flanders were selected.
}
\usage{data("HAV_BE_02")}
\format{
  A data frame with 5887 observations on the following 3 variables.
  \describe{
    \item{\code{birth.cohort}}{A numeric vector containing birth years of individuals.}
    \item{\code{age}}{A numeric vector containing individual ages.}
    \item{\code{status}}{A numeric vector containing individual HAV immunological status (boolean).}
  }
}
\details{
See page 62 of the book.
}
\examples{
data(HAV_BE_02)
}
\keyword{datasets}
\seealso{\code{\link{HAV_BE_9394}}}
