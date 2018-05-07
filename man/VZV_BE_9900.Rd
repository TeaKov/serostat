\name{VZV_BE_9900}
\alias{VZV_BE_9900}
\docType{data}
\title{
Varicella-Zoster Virus Serological Survey Data (Flanders, Belgium, 1999-2000)
}
\description{
The Varicella-Zoster virus, also known as human herpes virus 3, is one of eight herpes viruses known to affect humans. Primary VZV infection results in chickenpox (varicella). The age-specific seroprevalence of VZV antibodies was assessed in Flanders (Belgium) between October 1999 and April 2000. Sera from 1673 individuals, aged 1-44 years, were obtained.
}
\usage{data("VZV_BE_9900")}
\format{
  A data frame with 44 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{A numeric vector containing an ID number.}
    \item{\code{age}}{A numeric vector containing age of an individual.}
    \item{\code{pos}}{A numeric vector containing the number of seropositives in each age group.}
    \item{\code{neg}}{A numeric vector containing the number of seronegatives in each age group.}
    \item{\code{ntot}}{A numeric vector containing the total number of individuals in each age group.}
  }
}
\details{
See page 70 in the book.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/pubmed/12424583}
}
\examples{
data(VZV_BE_9900)
}
\keyword{datasets}
\seealso{\code{\link{VZV_BE_9900}}, \code{\link{VZV_B19_BE_0103}}}
