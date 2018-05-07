\name{HAV_BE_9394}
\alias{HAV_BE_9394}
\docType{data}
\title{
Hepatitis A Serological Survey Data (Flanders, Belgium, 1993-1994)
}
\description{
The purpose of this study, conducted in 1993 and early 1994, was to obtain data on the prevalence of Hepatitis A in Flanders, Belgium, and to analyze the epidemiological pattern of HAV. The dataset contains aggregated serological test results on 3161 Belgian individuals together with their age (mostly given in one year age intervals), ranging from 0.5 to 85 years.
}
\usage{data("HAV_BE_9394")}
\format{
  A data frame with 86 observations on the following 5 variables.
  \describe{
    \item{\code{Id}}{A numeric vector containing an ID number.}
    \item{\code{Age}}{A numeric vector containing age.}
    \item{\code{Pos}}{A numeric vector containing the number of seropositives in each age group.}
    \item{\code{Neg}}{A numeric vector containing the number of seronegatives in each age group.}
    \item{\code{Tot}}{A numeric vector containing the total number of individuals in each age group.}
  }
}
\details{
See page 62 of the book.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/pubmed/9258525}
}
\examples{
data(HAV_BE_9394)
}
\keyword{datasets}
\seealso{\code{\link{HAV_BE_02}}}
