\name{GLMPolyModels}
\alias{GLMPolyModels}
\alias{GLMPolyModels.default}
\alias{Muench}
\alias{MuenchAlt}
\alias{Griffiths}
\alias{GrenfellAnderson}
\alias{Weibull}
\title{Modeling the Force of Infection Using Polynomial Models}
\usage{
GLMPolyModels(parameters, fun)
}
\arguments{
\item{parameters}{Contains all input parameters: \bold{Age} describes the population's age categories, \bold{Pos} represent the number of seropositive individuals by age category, \bold{Tot} are the total number of individuals in each specified age category, and \bold{Bin} represents the response variable (binary, 1=past infection, 0=otherwise). One either specifies both the \bold{Pos} and \bold{Tot} parameters or the \bold{Bin} parameter. For the Griffith's model it is required to specify a \bold{tau} parameter, as the interval \eqn{(0,\tau)} represents the period of inherited immunity and since more than 95\% of the measles cases occur before the age of 10 years. The order of the polynomial of Grenfell's and Anderson's model can be specified by the \bold{order} parameter (default=2).}
\item{fun}{This specifies which polynomial model you want to use: \bold{Muench} is the catalytic model of Muench, \bold{MuenchAlt} is the alternative model of Muench, \bold{Griffiths} is Griffiths' model in which the force of infection increases linearly by age, and \bold{GrenfellAnderson} is Grenfell's and Anderson's model which extends the models of Muench and Griffiths to allow for higher-order polynomial functions representing force of infection. The \bold{Weibull} model implies that the force of infection is a monotone function of exposure time, and consequently the infection times are Weibull distributed.}
}
\description{
Models the (age-dependent) force of infection using polynomial generalized linear models, given input parameters and the specification of a model (see arguments). All models for the force of infection are cast into the well-known generalized linear model framework for binomial response data and are relying on different link-functions depending on the specified model.
}
\examples{
# Load Bulgarian Hepatitis A data
data("HAV_BUL_64")

# Prepare the parameters
params <- list(Tot=HAV_BUL_64$Tot, Pos=HAV_BUL_64$Pos, Age=HAV_BUL_64$Age)

# Muench's model
result <- GLMPolyModels(params, Muench)
summary(result$res)

# Muench's alternative model
result <- GLMPolyModels(params, MuenchAlt)
summary(result$res)

# Griffiths' model
params_griffiths <- c(params, tau=1)
result <- GLMPolyModels(params_griffiths, Griffiths)
summary(result$res)

# Grenfell's and Anderson's model
params_ga <- c(params, order=2)
result <- GLMPolyModels(params_ga, GrenfellAnderson)
summary(result$res)

# Load Belgian Hepatitis C data
data("HCV_BE_06_EXTRA")
data("HCV_BE_06")

# Prepare the parameters
params <- list(duration=HCV_BE_06_EXTRA$duration,
  infected=HCV_BE_06_EXTRA$infected, Tot=HCV_BE_06$V2,
  Pos=(HCV_BE_06$V2*HCV_BE_06$V4), Age=HCV_BE_06$V1)
result <- GLMPolyModels(params, Weibull)
summary(result$res)
}
