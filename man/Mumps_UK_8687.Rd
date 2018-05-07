\name{Mumps_UK_8687}
\alias{Mumps_UK_8687}
\docType{data}
\title{
Mumps Serological Survey Data (UK, 1986-1987)
}
\description{
Mumps is a childhood disease that occurs worldwide and is caused by a paramyxovirus. The most common symptoms of mumps are bilateral parotid swelling, fever, headache, and orchitis. This dataset contains a total of 8716 serum samples collected between November 1986 and December 1987 from five public health laboratories in different parts of the UK. The survey, covering subjects from 1 to over 65 years of age, provides information on the seroprevalence of mumps by age.
}
\usage{data("Mumps_UK_8687")}
\format{
  A data frame with 44 observations on the following 4 variables.
  \describe{
    \item{\code{age}}{A numeric vector containing age.}
    \item{\code{pos}}{A numeric vector containing the number of seropositives in each age group.}
    \item{\code{neg}}{A numeric vector containing the number of seronegatives in each age group.}
    \item{\code{ntot}}{A numeric vector containing the total number of individuals in each age group.}
  }
}
\details{
See page 66 of the book.
}
\examples{
data(Mumps_UK_8687)
}
\keyword{datasets}
\seealso{\code{\link{Rubella_UK_8687}}, \code{\link{Rubella_Mumps_UK_8687}}}
