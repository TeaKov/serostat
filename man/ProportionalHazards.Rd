\name{Proportional Hazards Model}
\alias{gam}
\alias{PropHazard}
\title{Estimating Age-Time Dependent Prevalence and Force of Infection from Serial Prevalence Data}
\usage{
gam(formula, family, \ldots)
}
\arguments{
\item{formula}{A GAM formula, or a list of formulae (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). These are exactly like the formula for a GLM except that smooth terms, \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} and \code{\link{t2}}, can be added to the right hand side to specify that the linear predictor depends on smooth functions of predictors (or linear functionals of these).}
\item{family}{This is a family object specifying the distribution and link to use in fitting etc (see \code{\link{glm}} and \code{\link{family}}). See \code{\link{family.mgcv}} for a full list of what is available, which goes well beyond exponential family. Note that quasi families actually result in the use of extended quasi-likelihood if method is set to a RE/ML method (McCullagh and Nelder, 1989, 9.6).}
\item{\ldots}{See \code{\link{gam}}.}
}
\description{
The use of serological surveys is nowadays a common way to study the epidemiology of many infections. In case a single cross-sectional survey is available, one needs to assume that the disease is in steady state. While reasonable for some infections, this assumption might be untenable for other situations. Models such as the proportional hazards model of Nagelkerke et al. (1999) estimate age- and time-specific prevalence and force of infection from a series of prevalence surveys.}
\examples{
library(mgcv)
# Load Dutch Tuberculosis Data
data("TB_NED_6673")
a <- TB_NED_6673$AGE
birthyear <- TB_NED_6673$BRTHYR
x <- birthyear-min(birthyear)
xf <- as.factor(birthyear)
g <- TB_NED_6673$SEX
ts <- a+birthyear
s <- TB_NED_6673$PPD
p <- s/TB_NED_6673$N
f <- TB_NED_6673$N-s
y <- cbind(s,f)

# Model 1 (parametric model) (Ades and Noke)
gamfit <- gam(y~a+x,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 2 (parametric model) (Ades and Noke)
gamfit=gam(y~a+x+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 3 (parametric model) (Ades and Noke)
gamfit=gam(y~a*x+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),sum(gamfit$edf),sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 4 (semi-parametric model) (Nagelkerke)
gamfit=gam(y~s(a)+x,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 5 (semi-parametric model) (Nagelkerke)
gamfit=gam(y~s(a)+x+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 6 (semi-parametric model) (Nagelkerke)
gamfit=gam(y~s(a)+x*g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 7 (non-parametric model)
gamfit=gam(y~s(a)+s(x),family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 8 (non-parametric model)
gamfit=gam(y~s(a)+s(x)+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 9 (non-parametric model)
gamfit=gam(y~s(a)+s(x)+s(x,by=g)+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 10 (non-proportional hazard)
gamfit=gam(y~s(a,x),family=binomial(link="cloglog"))
c(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 11 (non-proportional hazard)
gamfit=gam(y~s(a,x)+g,family=binomial(link="cloglog"))
cbind(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))

# Model 12 (non-proportional hazard)
gamfit=gam(y~s(a,x)+g+s(a,x,by=g),family=binomial(link="cloglog"))
c(sum(gamfit$edf),AIC(gamfit),AIC(gamfit, k = log(nrow(y))))
}
