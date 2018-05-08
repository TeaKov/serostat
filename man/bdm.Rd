\name{Bivariate Dale Model}
\alias{BDM}
\alias{vgam}
\title{Modeling the Prevalence and the Force of Infection Directly from Antibody Levels}
\usage{
vgam(formula, family, constraints, \ldots)
}
\arguments{
\item{formula}{A symbolic description of the model to be fit. The RHS of the formula is applied to each linear/additive predictor, and usually includes at least one s term. Different variables in each linear/additive predictor can be chosen by specifying constraint matrices.}
\item{family}{A symbolic description of the model to be fit. The RHS of the formula is applied to each linear predictor. Different variables in each linear predictor can be chosen by specifying constraint matrices.}
\item{constraints}{an optional list of constraint matrices. The components of the list must be named with the term it corresponds to (and it must match in character format exactly). There are two types of input: \dQuote{lm}-type and \dQuote{vlm}-type. The former is a subset of the latter. The former has a matrix for each term of the LM matrix. The latter has a matrix for each column of the VLM matrix. After fitting, the constraints extractor function may be applied; it returns the \dQuote{vlm}-type list of constraint matrices by default. If \dQuote{lm}-type are returned by constraints then these can be fed into this argument and it should give the same model as before. Each constraint matrix must have M rows, and be of full-column rank. By default, constraint matrices are the M by M identity matrix unless arguments in the family function itself override these values, e.g., parallel (see CommonVGAMffArguments). If constraints is used it must contain all the terms; an incomplete list is not accepted.}
\item{\ldots}{See \code{\link{vgam}}.}
}
\description{
A first marginal model that can be considered is the bivariate Dale model (BDM). The BDM relates the probability of past or current infection for both diseases to the age at infection.
}
\examples{
library(VGAM)
# Load Belgian B19 data
data("VZV_B19_BE_0103")
VZV_B19_BE_0103 <- VZV_B19_BE_0103[!is.na(VZV_B19_BE_0103$VZVres)&
  !is.na(VZV_B19_BE_0103$parvores)&!is.na(VZV_B19_BE_0103$age)&
  VZV_B19_BE_0103$age<70&VZV_B19_BE_0103$age>=1,]
VZV_B19_BE_0103 <- VZV_B19_BE_0103[order(VZV_B19_BE_0103$age),]

y1 <- VZV_B19_BE_0103$VZVres
y2 <- VZV_B19_BE_0103$parvores
age <- VZV_B19_BE_0103$age

a <- unique(age)
covariate <- seq(min(age),max(age),1)

# Counts per age-value
PP <- as.vector(hist(age[y1=="1"&y2=="1"], plot=FALSE, breaks=c(0,a))$counts)
PN <- as.vector(hist(age[y1=="1"&y2=="0"], plot=FALSE, breaks=c(0,a))$counts)
NP <- as.vector(hist(age[y1=="0"&y2=="1"], plot=FALSE, breaks=c(0,a))$counts)
NN <- as.vector(hist(age[y1=="0"&y2=="0"], plot=FALSE, breaks=c(0,a))$counts)

data <- data.frame(NN,NP,PN,PP,a=sort(a))

# BDM without constraints
\donttest{sel <- NULL
k.vec <- seq(0,0.1,0.001)
for (j in 1:length(k.vec)) {
  k <- k.vec[j]
  fit.s <- vgam(cbind(NN,NP,PN,PP)~s(a,spar=k),binom2.or(zero=NULL),data)
  BIC <- deviance(fit.s)+log(dim(data)[1])*(3*dim(data)[1]-df.residual(fit.s))
  sel <- rbind(sel,c(k,BIC))
}
sparopt<-sel[which.min(sel[,2]),1]

## Optimal value
fit.s1 <- vgam(cbind(NN,NP,PN,PP)~s(a,spar=sparopt),binom2.or(zero=NULL),data)
summary(fit.s1)
BIC1<-deviance(fit.s1)+log(dim(data)[1])*(3*dim(data)[1]-df.residual(fit.s1))


# Age-independent OR
fit.s2 <- vgam(cbind(NN,NP,PN,PP)~s(a,spar=sparopt),binom2.or(zero=NULL),
  constraints=list("(Intercept)"=diag(3),"s(a, spar = sparopt)"=diag(3)[,1:2]))
summary(fit.s2)
BIC2 <- deviance(fit.s2)+log(dim(data)[1])*(3*dim(data)[1]-df.residual(fit.s2))

# BDM Bootstrapping CI
runs <- 1000
runs.mat.pi1 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi2 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.or <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi00 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi01 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi10 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi11 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi2condpi1 <- matrix(NA,nrow=runs,ncol=length(covariate))
runs.mat.pi1condpi2 <- matrix(NA,nrow=runs,ncol=length(covariate))

for (i in 1:runs) {
  print(c("bootstrap",i))
  bootsample <- sample(c(1:length(a)),length(a),replace=TRUE)
  fit.boot <- vgam(cbind(NN,NP,PN,PP)~s(a,spar=sparopt),
    binom2.or(zero=NULL),data=data[bootsample,])
  runs.mat.pi1[i,] <-approx(x=a[bootsample],y=predictors(fit.boot)[,1],xout=covariate)$y
  runs.mat.pi2[i,] <-approx(x=a[bootsample],y=predictors(fit.boot)[,2],xout=covariate)$y
  runs.mat.or[i,] <-approx(x=a[bootsample],y=predictors(fit.boot)[,3],xout=covariate)$y
  runs.mat.pi00[i,] <-approx(x=a[bootsample],y=fitted(fit.boot)[,1],xout=covariate)$y
  runs.mat.pi01[i,] <-approx(x=a[bootsample],y=fitted(fit.boot)[,2],xout=covariate)$y
  runs.mat.pi10[i,] <-approx(x=a[bootsample],y=fitted(fit.boot)[,3],xout=covariate)$y
  runs.mat.pi11[i,] <-approx(x=a[bootsample],y=fitted(fit.boot)[,4],xout=covariate)$y
  runs.mat.pi1condpi2[i,]<-approx(x=a[bootsample],y=fitted(fit.boot)[,4]/
    (fitted(fit.boot)[,2]+fitted(fit.boot)[,4]),xout=covariate)$y
  runs.mat.pi2condpi1[i,]<-approx(x=a[bootsample],y=fitted(fit.boot)[,4]/
    (fitted(fit.boot)[,3]+fitted(fit.boot)[,4]),xout=covariate)$y
}}
}
\details{
See page 186 of the book.
}
\seealso{\code{\link{vgam}}}
