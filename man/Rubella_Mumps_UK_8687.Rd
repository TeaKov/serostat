\name{Rubella_Mumps_UK_8687}
\alias{Rubella_Mumps_UK_8687}
\docType{data}
\title{
Rubella and Mumps Serological Survey Data (UK, 1986-1987)
}
\description{
For feasibility and economical reasons, serum samples are often tested for more that one antigen. In this way, the (past) disease status of individuals on multiple diseases is known, and this allows studying the association in acquisition of several infections. Data from mumps and rubella are obtained from a large survey on the prevalence of antibodies to these infections in the UK, as described in \link{Rubella_UK_8687} and \link{Mumps_UK_8687}. In this survey, in total 8179 individuals were tested for both rubella and mumps IgG antibodies.
}
\usage{data("Rubella_Mumps_UK_8687")}
\format{
  A data frame with 44 observations on the following 5 variables.
  \describe{
    \item{\code{Age}}{A numeric vector containing age.}
    \item{\code{PP}}{A numeric vector containing the number of seropositives for both rubella and mumps.}
    \item{\code{PN}}{A numeric vector containing the number of individuals that are seropositive for rubella and seronegative for mumps.}
    \item{\code{NP}}{A numeric vector containing the number of individuals that are seropositive for mumps and seronegative for rubella.}
    \item{\code{NN}}{A numeric vector containing the number of seronegatives for both rubella and mumps.}
  }
}
\details{
See page 71 in the book.
}
\examples{
data(Rubella_Mumps_UK_8687)
}
\keyword{datasets}
\seealso{\code{\link{Rubella_UK_8687}}, \code{\link{Mumps_UK_8687}}}
