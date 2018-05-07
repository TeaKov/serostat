\name{Fitting WAIFW Structures - waifw.fitter}
\alias{waifw.fitter}
\title{Fitting WAIFW Structures}
\usage{
waifw.fitter(a, y, n, waifw.choice, muy, breaks, N, D, Lmax, A, startpar, plots)
}
\arguments{
\item{a}{Numeric vector of ages.}
\item{y}{The response variable (binary, 1=past infection, 0=otherwise).}
\item{n}{Number of individuals.}
\item{waifw.choice}{The choice of WAIFW structure.}
\item{muy}{Age-specific mortality rates.}
\item{breaks}{Breakpoints of the age-categories considered in the WAIFW approach.}
\item{N}{Population size.}
\item{D}{Mean duration of infectiousness.}
\item{Lmax}{Maximum age.}
\item{A}{Age of loss of maternal immunity (0<A<1).}
\item{startpar}{Starting values for the \dQuote{nlm} method used in this function.}
\item{plots}{Generate plots during fitting process (TRUE or FALSE.}
}
\description{
This function can be used to fit irregular WAIFW structures under constrained optimization. It requires the choice of WAIFW matrix as input (\dQuote{waifw.choice}) and uses as starting values the parameters from \code{\link{waifw.6parms}}.
}
\examples{
ND <- c(489,47,29,21,12,12,16,15,15,6,6,14,17,19,17,23,34,33,62,71,68,68,78,71,71,96,86,83,79,
  80,83,93,126,120,121,132,135,176,161,193,196,218,257,277,331,376,356,435,460,453,535,545,
  576,668,692,759,722,819,939,1015,1051,973,1113,996,940,1074,1252,1367,1468,1541,1661,1838,
  2012,2236,2517,2793,2938,2994,3311,3516,3727,3857,4088,4161,4261,4274,4061,2509,2049,2159,
  2205,2550,2330,1992,1569,1242,1000,726,533,996)
PS <- c(118366,117271,114562,113894,116275,118030,116761,117742,119583,119887,118963,119958,
  124637,129143,131030,129724,127187,126433,124377,124883,122201,124482,126459,130129,133897,
  135009,134516,133495,132705,132040,130602,135638,140537,146151,150467,152113,151656,151412,
  153371,158268,162456,167652,164871,161671,162060,159735,160672,157030,153820,151114,148978,
  145929,142374,141215,135525,135968,134692,135991,134291,134131,113024,112198,105880,92772,
  84462,93787,100820,101866,97208,94145,92451,93027,91640,93593,91933,89900,81718,77891,73104,
  70082,67057,62178,57642,51786,47466,42065,28004,17186,14492,13838,13957,13358,10442,8063,5604,
  4289,2843,2068,1368,2146)
AGE<-c(0:(length(ND)-1))

estimL <- estimateLifeExpectancy(ND, PS, AGE)

# Using Belgian B19 data
data("VZV_B19_BE_0103")
VZV_B19_BE_0103 <- VZV_B19_BE_0103[!is.na(VZV_B19_BE_0103$parvores)&
  !is.na(VZV_B19_BE_0103$age)&VZV_B19_BE_0103$age<70&VZV_B19_BE_0103$age>=1,]
VZV_B19_BE_0103 <- VZV_B19_BE_0103[order(VZV_B19_BE_0103$age),]

y <- VZV_B19_BE_0103$parvores
a <- VZV_B19_BE_0103$age

# Age category
breakpoints <- c(0.5,2,6,12,19,31,100)

# Mean duration of infectiousness
D <- 6/365

# Maximum life (if type mortality this is the life expectancy)
Lmax<-100

# Age of loss of maternal immunity (0<A<1)
A <- 0.5

# Mortality function
My <- estimL$My[1:Lmax]
muy <- estimL$muy[1:Lmax]

# Population size
N <- 10511382

# Fit WAIFW matrices W2 and W3
w2 <- waifw.fitter(a=a,y=y,waifw.choice=2,muy=muy,breaks=breakpoints,
  N=N,D=D,Lmax=Lmax,A=A,plots="TRUE")
w3 <- waifw.fitter(a=a,y=y,waifw.choice=3,muy=muy,breaks=breakpoints,
  N=N,D=D,Lmax=Lmax,A=A,plots="TRUE")
}
\seealso{\code{\link{pwcrate}}, \code{\link{waifw.6parms}}}
