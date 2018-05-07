\name{HCV_BE_06_EXTRA}
\alias{HCV_BE_06_EXTRA}
\docType{data}
\title{
Hepatitis C Serological Survey Data in IDUs (extra) (Belgium, 2006)
}
\description{
The Hepatitis C virus (HCV) is the leading cause of known liver disease in most industrialized countries. Mathei et al. (2006) presented a study of HCV infection among injected drug users (IDUs). The data consists of IDUs from three areas in Belgium: Charleroi, Antwerp, and the province of Limburg.
}
\usage{data("HCV_BE_06_EXTRA")}
\format{
  A data frame with 421 observations on the following 2 variables.
  \describe{
    \item{\code{duration}}{A numeric vector containing the duration of injection.}
    \item{\code{infected}}{A numeric vector containing the immunological status with regard to HCV infection (boolean).}
  }
}
\details{
See pages 65 and 97 of the book.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/pubmed/16901287}
}
\examples{
data(HCV_BE_06_EXTRA)
}
\keyword{datasets}
