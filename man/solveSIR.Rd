\name{solveSIR}
\alias{solveSIR}
\alias{solveSIR.default}
\alias{SIR_NoVac}
\alias{SIR_Vac}
\alias{SIR_ConstFOI}
\alias{SIR_SubPop}
\title{Solving SIR models}
\usage{
solveSIR(t, parameters, state, fun)
}
\arguments{
\item{t}{Timings of the integration.}
\item{parameters}{Contains all the model parameters: natural mortality rate \bold{mu}, disease-related mortality \bold{alpha}, transmission parameter \bold{beta}, recovery rate \bold{nu}, portion of vaccinated newborns \bold{p}, and total population size \bold{N}. The parameter \bold{ndims} specifies the number of row and columns of the subpopulation matrix, and \bold{nparams} specifies the dimension of the `nu' parameter (SIR_SubPop).}
\item{state}{The initial values of the model: \bold{S} are the susceptibles, \bold{I} are the infected, and \bold{R} are the immune or recovered.}
\item{fun}{This specifies which SIR model you want to use: \bold{SIR_Vac} is the SIR model with a portion of vaccinated newborns (Halloran), \bold{SIR_NoVac} is the SIR model without any vaccinated portion, \bold{SIR_ConstFOI} is the SIR model with constant force of infection at endemic state and without a portion of vaccinated newborns, and \bold{SIR_SubPop} is the SIR model that uses a mixing matrix that specifies a subpopulation.}
}
\description{
Solves a mathematical SIR model, with or without vaccinated portion of newborns or constant force of infection, given timings, parameters and state.
}
\examples{
# Without portion of vaccinated newborns
t <- seq(0,1600,0.01)
state <- c(S=4999,I=1,R=0)
params <- c(mu=1/75,beta=0.0005,nu=1,N=5000,alpha=0)
result <- solveSIR(t, params, state, SIR_NoVac)
print(result)

# With portion of vaccinated newborns
t <- seq(0,1600,0.01)
state <- c(S=4999,I=1,R=0)
params <- c(mu=1/75,beta=0.0005,nu=1,N=5000,p=0.4,alpha=0)
result <- solveSIR(t, params, state, SIR_Vac)
print(result)

# With constant force of infection
t <- seq(0,90,by=0.01)
state <- c(s=0.99,i=0.01,r=0)
params <- c(lambda = 0.05, nu=1/(14/365))
result <- solveSIR(t, params, state, SIR_ConstFOI)
print(result)

# With subpopulation mixing matrix
t <- seq(0,10000,by=0.01)
state <- c(s=c(0.8,0.8),i=c(0.2,0.2), r=c(0.0,0.0))
params <- c(betatilde=c(0.05,0.0,0.0,0.05), nu=c(1/30,1/30), mu=0.001, ndims=2, nparams=2)
result <- solveSIR(t, params, state, SIR_SubPop)
print(result)
}
