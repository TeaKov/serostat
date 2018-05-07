\name{NonLinearModels}
\alias{NonLinearModels}
\alias{NonLinearModels.default}
\alias{FarringtonABG}
\alias{FarringtonAB_G}
\title{Modeling the Force of Infection using Non-Linear Models}
\usage{
NonLinearModels(data, parameters, fun)
}
\arguments{
\item{data}{Contains all input parameters: \bold{Age} describes the population's age categories, \bold{Pos} represents the number of seropositive individuals by age category, and \bold{Tot} are the total number of individuals in each specified age category.}
\item{parameters}{Specifies the parameters used in the non-linear force of infection equation: \eqn{\lambda(a) = (a\alpha - \gamma)e^{-a\beta} + \gamma}, where \bold{gamma} is the long-term residual value of the force of infection (Farrington et al., 1990). In order to ensure that the force of infection satisfies \eqn{\lambda(a) \ge 0}, the parameter space is constrained to be non-negative.}
\item{fun}{This specifies which non-linear model for the force of infection you want to use: \bold{FarringtonABG} is Farrington's model where alpha, beta, and gamma are used as model parameters, \bold{FarringtonAB_G} is an alternative simplified model in which gamma is kept as a fixed parameter and alpha and beta are the model parameters that are estimated.}
}
\description{
Models the (age-dependent) force of infection using non-linear models, given input parameters and the specification of a model (see arguments).
}
\examples{
# Load UK's Rubella data from '86-'87
data("Rubella_UK_8687")

# Prepare the data
data <- list(Tot=(Rubella_UK_8687$Pos+Rubella_UK_8687$Neg),
  Pos=Rubella_UK_8687$Pos, Age=Rubella_UK_8687$Age)

# Prepare the parameters
params <- list(alpha=0.07, beta=0.1, gamma=0.03)

# Farrington's model, with alpha, beta, and gamma as starting values for MLE
result <- NonLinearModels(data, params, FarringtonABG)
summary(result)

# Farrington's model, with alpha and beta as starting values for MLE, and gamma as a fixed parameter
params <- list(alpha=0.07, beta=0.1, gamma=0)
result <- NonLinearModels(data, params, FarringtonAB_G)
summary(result)
}
\seealso{\cite{Farrington, C. P. (1990). Modelling forces of infection for measles, mumps and rubella, 9(8), 953-967.}}
