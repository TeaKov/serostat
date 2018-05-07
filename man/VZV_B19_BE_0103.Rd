\name{VZV_B19_BE_0103}
\alias{VZV_B19_BE_0103}
\docType{data}
\title{
Varicella-Zoster Virus and Parvovirus B19 Serological Survey Data (Belgium, 2001-2003)
}
\description{
For feasibility and economical reasons, serum samples are often tested for more that one antigen. In this way, the (past) disease status of individuals on multiple diseases is known, and this allows studying the association in acquisition between several infections. Serological data on Varicella-Zoster Virus (VZV) and Parvovirus B19 (B19) infections are collected in a serosurvey in Belgium, as described in \code{\link{VZV_BE_9900}} and \code{\link{B19_ALL}}. In total, 3080 sera were tested for VZV and 2657 sera were tested for B19, from which 2382 sera were tested for both VZV and Parvovirus B19. This dataset can be used to, for example, visualize individuals that tested seropositive for both VZV and B19, those that only tested seropositive for VZV, those that only tested seropositive for B19, and those that tested seronegative for both.
}
\usage{data("VZV_B19_BE_0103")}
\format{
  A data frame with 3374 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{A factor encoding a unique identification number for each individual.}
    \item{\code{sex}}{A numeric vector containing gender (female: 1, male: 2).}
    \item{\code{age}}{A numeric vector containing age of an individual.}
    \item{\code{parvouml}}{A numeric vector containing the continuous IgG antibody concentration against Parvovirus B19 (in ml).}
    \item{\code{parvores}}{A numeric vector containing the immunological status with regard to Parvovirus B19.}
    \item{\code{VZVmUIml}}{A numeric vector containing the continuous IgG antibody concentration against VZV (in ml).}
    \item{\code{VZVres}}{A numeric vector containing the immunological status with regard to VZV.}
  }
}
\details{
See page 71 in the book.
}
\examples{
data(VZV_B19_BE_0103)
}
\keyword{datasets}
\seealso{\code{\link{VZV_BE_9900}}, \code{\link{B19_ALL}}}
