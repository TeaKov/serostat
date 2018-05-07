\name{Semiparametric Regression}
\alias{SemiParametricRegression}
\alias{SemiParametricRegression.default}
\alias{SmoothSplines}
\alias{CubicRegSplines}
\alias{CubicRegSplinesMM}
\alias{ThinPlateRegSplines}
\alias{ThinPlateRegSplinesMM}
\alias{AdaptiveSplineSmoothing}
\title{Semiparametric Approaches to Model the Prevalence and Force of Infection}
\usage{
SemiParametricRegression(parameters, fun)
}
\arguments{
\item{parameters}{Contains all the parameters of a model: \bold{y} represents a vector of the immunological statuses of all individuals with regard to the pathogen under consideration, \bold{a} is a vector of the ages of all individuals at the time of data collection, and \bold{s} is a vector containing information about the individual's sex (as categorical variable). On top of that, approach-specific parameters should be specified. \bold{SmoothSplines} requires the specification of the \bold{dfopts} parameter (the target equivalent degrees of freedom, used as a smoothing parameter) and \bold{link} parameter (``logit'' or ``cloglog''). \bold{CubicRegSplines}, \bold{ThinPlateRegSplines}, \bold{CubicRegSplinesMM}, and \bold{ThinPlateRegSplinesMM} require that one specifies the \bold{link} parameter (``logit'' or ``cloglog'').}
\item{fun}{This specifies which semiparametric regression model you want to use: \bold{SmoothSplines} avoids knot selection problems by using the maximal set of knots. \bold{CubicRegSplines} joins (cubic) polynomials at the knots of the spline to ensure continuity and differentiablility up to degree two. \bold{ThinPlateRegSplines} does not use knots, is computationaly harder but provides nested models which for building models is in line with general linear modeling methods. \bold{CubicRegSplinesMM} and \bold{ThinPlateRegSplinesMM} are the same as \bold{CubicRegSplines} and \bold{ThinPlateRegSplines}, but model parameter estimation is performed using the GLMM framework and the PQL-estimation method. \bold{AdaptiveSplineSmoothing} employs the Laplace approximation to fit an adaptive penalized spline in the GLMM framework.}
}
\description{
Models the force of infection using semi-parametric regression models, given input parameters and the specification of a model (see arguments). For more details on the B-splines approach, go to help("pspline.fit").
}
\examples{
# Load Belgian B19 data
data("VZV_B19_BE_0103")
subset<-(VZV_B19_BE_0103$age>0.5)&(VZV_B19_BE_0103$age<76)&(!is.na(VZV_B19_BE_0103$age))&
  !is.na(VZV_B19_BE_0103$parvores)
VZV_B19_BE_0103<-VZV_B19_BE_0103[subset,]
y<-VZV_B19_BE_0103$parvores[order(VZV_B19_BE_0103$age)]
a<-VZV_B19_BE_0103$age[order(VZV_B19_BE_0103$age)]

# Smoothing splines, with logit link-function
dfopts <- 3.5
params <- list(y=y, a=a, dfopts=dfopts, link="logit")
result <- SemiParametricRegression(params, SmoothSplines)

# Smoothing splines, with cloglog link-function
dfopts <- 4.5
params <- list(y=y, a=a, dfopts=dfopts, link="cloglog")
result <- SemiParametricRegression(params, SmoothSplines)

# Cubic regression splines, with logit link-function
params <- list(y=y, a=a, link="logit")
result <- SemiParametricRegression(params, CubicRegSplines)

# Cubic regression splines, with cloglog link-function
params <- list(y=y, a=a, link="cloglog")
result <- SemiParametricRegression(params, CubicRegSplines)

# Thin Plate regression splines, with logit link-function
params <- list(y=y, a=a, link="logit")
result <- SemiParametricRegression(params, ThinPlateRegSplines)

# Thin Plate regression splines, with cloglog link-function
params <- list(y=y, a=a, link="cloglog")
result <- SemiParametricRegression(params, ThinPlateRegSplines)

# Cubic regression splines, with logit link-function and gamm
params <- list(y=y, a=a, link="logit")
result <- SemiParametricRegression(params, CubicRegSplinesMM)

# Cubic regression splines, with cloglog link-function and gamm
params <- list(y=y, a=a, link="cloglog")
result <- SemiParametricRegression(params, CubicRegSplinesMM)

# Thin Plate regression splines, with logit link-function and gamm
params <- list(y=y, a=a, link="logit")
result <- SemiParametricRegression(params, ThinPlateRegSplinesMM)

# Thin Plate regression splines, with cloglog link-function and gamm
params <- list(y=y, a=a, link="cloglog")
result <- SemiParametricRegression(params, ThinPlateRegSplinesMM)

# Adaptive Spline Smoothing
library(AdaptFit)
kn.mean <- default.knots(a, 20)
kn.var <- default.knots(a,5)
params <- list(y=y, a=a, adap=TRUE, knots=kn.mean, var.knot=kn.var)
result <- SemiParametricRegression(params, AdaptiveSplineSmoothing)

# Adaptive Spline Smoothing
library(AdaptFit)
kn.mean <- default.knots(a, 20)
kn.var <- default.knots(a,5)
params <- list(y=y, a=a, adap=FALSE, knots=kn.mean)
result <- SemiParametricRegression(params, AdaptiveSplineSmoothing)
}
