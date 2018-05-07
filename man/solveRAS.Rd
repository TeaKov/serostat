\name{solveRAS}
\alias{solveRAS}
\alias{solveRAS.default}
\alias{RAS_FOI}
\alias{RAS_BETAS}
\alias{RAS_PDE}
\title{Solving Realistic Age-Structured SIR Models}
\usage{
solveRAS(config, parameters, state, fun, \ldots)
}
\arguments{
\item{config}{Contains the configuration of the simulation: \bold{Tinit} contains the timings of integration, \bold{Truns} specifies the amount of simulation runs to be executed, \bold{tol} is negates negative values in \dQuote{states}, \bold{res} specifies the number of subdivisions with respect to time, and \bold{cohort.size} specifies the cohort size.}
\item{parameters}{Contains all the parameters of the model: \bold{n.ages} is number of age categories, \bold{mu} is natural-related deaths, \bold{beta} is the WAIFW or social contact matrix (transmission parameter), \bold{FOI} is the age-dependent force of infection (transmission parameter), and \bold{alpha.rate} is the recovery rate.}
\item{state}{The initial values of the model: \bold{S} are the number of susceptibles, \bold{I} are the number of infected individuals, and \bold{R} are the number of immune or recovered individuals.}
\item{fun}{This specifies which realistic age-structured SIR model you want to use: \bold{RAS_FOI} uses an age-dependent force of infection (meaning that the \dQuote{FOI} parameter should be specified and not the \dQuote{betas} parameter), \bold{RAS_BETAS} is the SIR model that uses a WAIFW matrix or a social contact matrix as transmission parameter (meaning that the \dQuote{betas} parameter should be specified and not the \dQuote{FOI} parameter),  and \bold{RAS_PDE} also uses a WAIFW or social contact matrix as transmission parameter, but the model itself is a PDE model instead of an ODE model. Note that also here one should only use the \dQuote{betas} parameter instead of the \dQuote{FOI} parameter.}
\item{\ldots}{See \code{\link{lsoda}}}
}
\description{
Solves realistic age-structured SIR models, given an age-dependent force of infection, a WAIFW or social contact matrix or the latter that is implemented using a system of PDEs.
}
\examples{
# Mortality function for BE anno 2006
# The number of deaths during the year 2006 and the average population anno 2006-2007
ND <- c(489,47,29,21,12,12,16,15,15,6,6,14,17,19,17,23,34,33,62,71,68,68,78,71,71,96,86,83,79,
  80,83,93,126,120,121,132,135,176,161,193,196,218,257,277,331,376,356,435,460,453,535,545,576,
  668,692,759,722,819,939,1015,1051,973,1113,996,940,1074,1252,1367,1468,1541,1661,1838,2012,2236,
  2517,2793,2938,2994,3311,3516,3727,3857,4088,4161,4261,4274,4061,2509,2049,2159,2205,2550,2330,
  1992,1569,1242,1000,726,533,996)
PS2006 <- c(118366,117271,114562,113894,116275,118030,116761,117742,119583,119887,118963,119958,
  124637,129143,131030,129724,127187,126433,124377,124883,122201,124482,126459,130129,133897,
  135009,134516,133495,132705,132040,130602,135638,140537,146151,150467,152113,151656,151412,
  153371,158268,162456,167652,164871,161671,162060,159735,160672,157030,153820,151114,148978,
  145929,142374,141215,135525,135968,134692,135991,134291,134131,113024,112198,105880,92772,
  84462,93787,100820,101866,97208,94145,92451,93027,91640,93593,91933,89900,81718,77891,73104,
  70082,67057,62178,57642,51786,47466,42065,28004,17186,14492,13838,13957,13358,10442,8063,5604,
  4289,2843,2068,1368,2146)
PS2007 <- c(121718,119795,118426,115497,114720,117067,118696,117411,118410,120276,120530,119564,
  120635,125230,129754,131590,130406,128061,127594,125749,126481,124131,126329,128238,131953,
  135668,136899,136289,135193,134314,133529,132009,136806,141763,147274,151465,153140,152332,
  152104,153956,158802,162872,168005,165246,161831,162098,159919,160630,156922,153601,150843,
  148660,145653,141968,140728,134950,135309,134021,135186,133339,133033,112049,111107,104865,
  91833,83476,92583,99546,100390,95774,92563,90691,91139,89531,91210,89332,86949,78757,74757,
  69725,66407,63301,58200,53457,47599,43154,37679,24848,14954,12433,11668,11580,10843,8262,6242,
  4243,3076,2039,1414,2240)
PS <- (PS2006+PS2007)/2

# Initial parameters
Tinit <- c(1, 200)
Truns <- 200
months.days <- c(30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31)
cohort.size <- 120000
year.days <- sum(months.days)       # Length of the year
alpha.days <- 7                     # Infectious period (days)
alpha.rate <- year.days/alpha.days  # Infectiousness rate (years^-1)
n.ageclass = 100                    # Number of age classes
mu <- -log(1-(ND/PS)[1:100])        # Mortality rate for Belgium
tol <- 1.00E-10                     # Precision
resolution <- 365                   # Resolution (nr subdivisions wrt time)
p <- 0.75                           # Vaccination coverage
q <- 0.173                          # Proportionality factor for the contact matrix option 1

config <- list(Tinit=Tinit, Truns=Truns, tol=tol, res=resolution, cohort.size=cohort.size)

# Constant FOI
# Read in the age-dependent force of infection
data("FOI_VZV")

S <- rep(0, 100)
I <- rep(0, 100)
R <- rep(0, 100)

S <- (exp(-cumsum(FOI_VZV)))*cohort.size*exp(-cumsum(mu))
I <- 1/alpha.rate*FOI_VZV*exp(-cumsum(FOI_VZV))*cohort.size*exp(-cumsum(mu))
R <- cohort.size*exp(-cumsum(mu))-S-I

states <- c(S=S, I=I, R=R)
params <- c(n.ages=n.ageclass, mu=mu, foi=FOI_VZV, alpha.rate=alpha.rate)

res_foi <- solveRAS(config, params, states, RAS_FOI)

# WAIFW matrix
beta1 <- 1.413; beta2 <- 1.335; beta3 <- 1.064; beta4 <- 0.0; beta5 <- 0.343; beta6 <- 0.0
W <- matrix(beta6,ncol=100,nrow=100)
W[1:31,1:31] <- beta5; W[1:19,1:19] <- beta4; W[1:12,1:12] <- beta3
W[1:6,1:6] <- beta1; W[3:6,3:6] <- beta2; W<-W*(10^-4)

S <- rep(0, 100)
I <- rep(0, 100)
R <- rep(0, 100)

S <- (exp(-cumsum(FOI_VZV)))*cohort.size*exp(-cumsum(mu))
I <- 1/alpha.rate*FOI_VZV*exp(-cumsum(FOI_VZV))*cohort.size*exp(-cumsum(mu))
R <- cohort.size*exp(-cumsum(mu))-S-I

states <- c(S=S, I=I, R=R)
params <- c(n.ages=n.ageclass, mu=mu, betas=W, alpha.rate=alpha.rate)

res_waifw <- solveRAS(config, params, states, RAS_BETAS)

# Social contact data
# Use social contact data "close+15"
data("scd_close_p15")

cij <- matrix(min(scd_close_p15),nrow=100,ncol=100)
cij[1:86,1:86] <- as.matrix(scd_close_p15)
betas <- 365*q*cij

S <- rep(0, 100)
I <- rep(0, 100)
R <- rep(0, 100)

S <- (exp(-cumsum(FOI_VZV)))*cohort.size*exp(-cumsum(mu))
I <- 1/alpha.rate*FOI_VZV*exp(-cumsum(FOI_VZV))*cohort.size*exp(-cumsum(mu))
R <- cohort.size*exp(-cumsum(mu))-S-I

states <- c(S=S, I=I, R=R)
params <- c(n.ages=n.ageclass, mu=mu, betas=betas, alpha.rate=alpha.rate)

res_contact <- solveRAS(config, params, states, RAS_BETAS)

# With PDEs
res_contact_pde <- solveRAS(config, params, states, RAS_PDE)
}
