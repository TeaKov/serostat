\name{B19_ALL}
\alias{B19_ALL}
\docType{data}
\title{
Parvovirus B19 Serological Survey Data (Belgium, England and Wales, Finland, Italy, Poland)
}
\description{
Parvovirus B19 is the infectious agent that causes erythema infectiosum, commonly known as slapped cheek syndrome or fifth disease. Seroprevalence surveys testing for parvovirus B19 IgG antibodies were performed on large representative national serum banks in Belgium, England and Wales, Finland, Italy, and Poland. The sera were collected between 1995 and 2004 and were obtained from residual sera submitted for routine laboratory testing.
}
\usage{data("B19_ALL")}
\format{
  A data frame with 13446 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{A factor encoding unique ID numbers for all individuals.}
    \item{\code{sex}}{A factor encoding the sex of each individual.}
    \item{\code{age}}{A numeric vector including the individual ages.}
    \item{\code{parvouml}}{A numeric vector including the individual volume of Parvovirus IgG antibodies (in ml).}
    \item{\code{parvores}}{A numeric vector (boolean) including the individual Parvovirus immunological responses based on threshold values pre-specified by the manufacturer of the diagnostic tests.}
    \item{\code{year}}{A factor encoding years of data collection.}
    \item{\code{country}}{A factor encoding countries (be: Belgium, ew: England and Wales, fi: Finland, it: Italy, pl: Poland).}
  }
}
\details{
See page 68 of the book.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870907/}
}
\examples{
data(B19_ALL)
}
\keyword{datasets}
\seealso{\code{\link{VZV_BE_9900}}, \code{\link{VZV_B19_BE_0103}}}
