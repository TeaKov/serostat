\name{HBV_Incidence_99}
\alias{HBV_Incidence_99}
\docType{data}
\title{
Hepatitis B Serological Survey Data (St.-Petersburg, Russia, 1999)
}
\description{
Data from a seroprevalence study conducted in St.-Petersburg (Russia) in 1999 (Mukomolov et al., 2000). The study intended to collect sera from 100 healthy persons (50 males, 50 females) in each of the following age groups: < 1, 1 - 2, 3 - 6, 7 - 10, 11 - 14, 15 - 19, 20 - 24, 25 - 29, 30 - 39, > 40 years of age (total original sample: 1003 sera). For the youngest age groups, the sample was taken primarily among children in the kindergartens and schools, and supplemented with sera from children entering hospitals with acute noninfectious pathologies (e.g. trauma, emergency surgery, pneumonia). Sera from teenagers and young adults were obtained from a variety of schools (excluding nursing schools) as well as from hospitalized persons requiring urgent surgery. Sera from adults were obtained from primary blood donors, pregnant women and persons entering the hospitals with acute noninfectious pathologies. All sera were tested for anti-HBc.
}
\usage{data("HBV_Incidence_99")}
\format{
  A data frame with 100 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{A numeric vector containing age.}
    \item{\code{men}}{A numeric vector containing the number of men per age category.}
    \item{\code{women}}{A numeric vector containing the number of women per age category.}
  }
}
\details{
See page 64 of the book.
}
\examples{
data(HBV_Incidence_99)
}
\keyword{datasets}
\seealso{\code{\link{HBV_Incidence_99_Dates}}, \code{\link{HBV_Population}}}
