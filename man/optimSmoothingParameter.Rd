\name{Penalized B-splines with Monotonicity Constraint - Finding Optimal Smoothing Parameter}
\alias{optimSmoothingParameter}
\title{Penalized B-splines with Monotonicity Constraint - Finding Optimal Smoothing Parameter}
\usage{
optimSmoothingParameter(mpspline.params, alpha.init, useAIC, upper)
}
\arguments{
\item{mpspline.params}{List of \code{mpspline.fit} function parameters. Required: the response variable \code{response}, and the explanatory variable on the abcissae \code{x.var}. See \code{\link{mpspline.fit}} for a full list of parameters.}
\item{alpha.init}{Initial value for the smoothness regularizing parameter (>= 0). Default=2.}
\item{useAIC}{Specifies if the AIC criteria is used for optimization. Otherwise, the BIC is used. Default=TRUE.}
\item{upper}{The upper bound on the variable that will be optimized by \code{\link{optim}}.}
}
\description{
This function tries to find the optimal smoothing parameter value used in \code{\link{mpspline.fit}}.
}
\examples{
# Using Belgian B19 data
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age<40.5)&(!is.na(VZV_B19_BE_0103$age))&
  (!is.na(VZV_B19_BE_0103$VZVmUIml))
data <- VZV_B19_BE_0103[subset,]

# Data to use when taking the continuous levels
z<-log(data$VZVmUIml[order(data$age)]+1)
a<-data$age[order(data$age)]

params <- list(response=z, x.var=a, ps.intervals=20, degree=3,
  link="identity", family="gaussian", order=2, kappa=1e8)

optimSmoothingParameter(params)

# Proportion positive, as a function of the corresponding half-year
# age categories, overlaid with the monotone p-spline fit and the FOI:
# with BIC optimal smoothing parameter
#
# Using Belgian B19 data
data("VZV_B19_BE_0103")
subset <- (VZV_B19_BE_0103$age<40.5)&(!is.na(VZV_B19_BE_0103$age))&
  (!is.na(VZV_B19_BE_0103$VZVmUIml)&(!is.na(VZV_B19_BE_0103$VZVres)))
data <- VZV_B19_BE_0103[subset,]

# Data to use when taking the binary indicators (different because some are inconclusive)
y <- data$VZVres[order(data$age)]
a <- data$age[order(data$age)]

params <- list(response=y, x.var=a, ps.intervals=20, degree=3,
  order=2, link="logit", family="binomial", kappa=1e8)

optimSmoothingParameter(params,useAIC=FALSE)
}
\seealso{\code{\link{mpspline.fit}}, \code{\link{optim}}}
