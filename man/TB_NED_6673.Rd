\name{TB_NED_6673}
\alias{TB_NED_6673}
\docType{data}
\title{
Tuberculosis Serological Survey Data (The Netherlands, 1966-1973)
}
\description{
Tuberculosis (TB) is a bacterial infection caused by Mycobacterium tuberculosis. TB can attack any organ, but most infections are restricted to the lungs. In 1966-1973, a study of tuberculosis was conducted in the Netherlands.
}
\usage{data("TB_NED_6673")}
\format{
  A data frame with 110 observations on the following 5 variables.
  \describe{
    \item{\code{PPD}}{A numeric vector containing the number of seropositives for tuberculosis.}
    \item{\code{N}}{A numeric vector containing the total number of individuals of that age category and sex.}
    \item{\code{AGE}}{A numeric vector containing age.}
    \item{\code{BRTHYR}}{A numeric vector containing the birth year (ranging 66-73).}
    \item{\code{SEX}}{A numeric vector containing sex (boolean, 0: female, 1: male).}
  }
}
\details{
See page 69 of the book.
}
\source{
\url{http://www.sciencedirect.com/science/article/pii/0041387983900211} and \url{http://www.ncbi.nlm.nih.gov/pubmed/10070676}
}
\examples{
data(TB_NED_6673)
}
\keyword{datasets}
