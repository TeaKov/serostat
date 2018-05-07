# This file contains polynomial models introduced in chapter 6

plot_sero_foi <- function(title, p, fitted.vals, coeffs, X) {
  par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=1,mgp=c(2, 0.5, 0),mar=c(4.1,4.1,4.1,4.1))
  plot(p$Age,p$Pos/p$Tot,cex=0.1*p$Tot,xlab="age",xlim=c(0,max(p$Age)),ylim=c(0,1),ylab="seroprevalence", main=paste("Seroprevalence and FOI,",title))
  lines(p$Age,1-fitted.vals,lwd=2)

  lines(p$Age,5*X%*%coeffs,lwd=2)
  axis(side=4,at=c(0.0,0.2,0.4),labels=c(0.00,0.04,0.08))
  mtext(side=4,"force of infection", las=3,line=2)
}

Muench <- function(p)
{
  res <- stats::glm(base::cbind(p$Tot-p$Pos, p$Pos)~-1+p$Age, family=stats::binomial(link="log"))

  X <- -matrix(rep(1,length(p$Age)))
  plot_sero_foi("Muench", p, res$fitted.values, res$coefficients, X)

  res <- list(X=X, res=res)
  res
}

MuenchAlt <- function(p)
{
  res <- stats::glm(base::cbind(p$Pos, p$Tot-p$Pos)~1, offset=base::log(p$Age), family=stats::binomial(link="cloglog"))

  X <- -matrix(rep(1,length(p$Age)))
  plot_sero_foi("MuenchAlt", p, res$fitted.values, res$coefficients, X)

  res <- list(X=X, res=res)
  res
}

Griffiths <- function(p)
{
  lb_g <- which(p$Age == p$tau)
  ub_g <- which(p$Age == 10)
  Tot_g <- p$Tot[lb_g:ub_g]
  Age_g <- p$Age[lb_g:ub_g]
  Pos_g <- p$Pos[lb_g:ub_g]

  res <- stats::glm(base::cbind(Tot_g-Pos_g, Pos_g)~-1+Age_g+base::I(Age_g^2), family=stats::binomial(link="log"))

  X <- -cbind(rep(1,length(Age_g)),2*Age_g)
  p_new <- list(Tot=Tot_g, Age=Age_g, Pos=Pos_g)
  plot_sero_foi("Griffiths", p_new, res$fitted.values, res$coefficients, X)

  res <- list(X=X, res=res)
  res
}

GrenfellAnderson <- function(p)
{
  Age_ga <- c("p$Age")
  X <- cbind(-rep(1,length(p$Age)))
  for(i in 2:(p$order+1)) {
    Age_ga <- c(Age_ga, sprintf("base::I(p$Age^%d)", i))
    X <- cbind(X, -i*(p$Age^(i-1)))
  }

  form_str <- paste("base::cbind(p$Tot-p$Pos, p$Pos)~-1+", paste(Age_ga, collapse="+"))

  res <- stats::glm(as.formula(form_str), family=stats::binomial(link="log"))
  plot_sero_foi("GrenfellAnderson", p, res$fitted.values, res$coefficients, X)

  res <- list(X=X, res)
  res
}

Weibull <- function(p)
{
  log.d <- base::log(p$duration)
  res <- stats::glm(p$infected~log.d, family=stats::binomial(link="cloglog"))

  plot(p$Age,p$Pos/p$Tot,cex=0.1*p$Tot,xlab="age",xlim=c(0,25),ylim=c(0,1),ylab="seroprevalence", main="Seroprevalence and FOI, Weibull")
  b0 <- stats::coef(res)[1]
  b1 <- stats::coef(res)[2]
  fitted <- stats::predict(res,type="response")
  fittedorig <- 1-base::exp(-base::exp(b0)*p$duration^b1)
  lines(sort(exp(log.d)),fittedorig[base::order(p$duration)],lwd=2)

  foi <- base::exp(b0)*b1*p$duration^(b1-1)
  foi.t <- 0.5*foi[p$duration>=0.5]
  t.t <- p$duration[p$duration>=0.5]
  lines(sort(t.t),foi.t[order(t.t)],lwd=2)
  axis(side=4,at=c(0.0,0.1,0.2),labels=c(0.0,0.2,0.4))
  mtext(side=4,"force of infection", las=3,line=1.5)

  res <- list(res=res)
  res
}

GLMPolyModels <- function(parameters, fun)
{
  #UseMethod("GLMPolyModels", fun)
  UseMethod("GLMPolyModels")
}

GLMPolyModels.default <- function(parameters, fun)
{
  function_string = base::as.character(substitute(fun))

  if(function_string == "Weibull")
  {
    noDuration <- !("duration" %in% base::names(parameters))
    noInfected <- !("infected" %in% base::names(parameters))

    if(noDuration)
      warning("d not specified!")
    if(noInfected)
      warning("infected not specified!")
    if(noDuration || noInfected)
      stop("Please specify the necessary parameters!")
  }
  else
  {
    noTot <- !("Tot" %in% base::names(parameters))
    noPos <- !("Pos" %in% base::names(parameters))
    noAge <- !("Age" %in% base::names(parameters))

    noBin <- !("Bin" %in% base::names(parameters))

    missingParameters <- FALSE

    if(noAge) {
      warning("Age not specified!")
      missingParameters <- TRUE
    }

    if(noBin && (noTot || noPos)) {
      if(noTot)
        warning("Tot not specified!")
      if(noPos)
        warning("Pos not specified!")

      missingParameters <- TRUE
    }

    if(noTot && noPos && noBin) {
      warning("Bin not specified!")
      missingParameters <- TRUE
    }

    if(function_string == "Griffiths" && !("tau" %in% base::names(parameters))) {
      warning("Tau not specified!")
      missingParameters <- TRUE
    }

    if(missingParameters)
      stop("Please specify the necessary parameters!")

    if(function_string == "Griffiths" && parameters$tau >= 10)
      stop("Tau should be less than 10!")

    if(function_string == "GrenfellAnderson") {
      if(!("order" %in% base::names(parameters))) {
        message("Assuming polynomial order 2...")
        parameters <- c(parameters, order=2)
      }
      else if(parameters$order < 1) {
        stop("Order should be larger or equal to 1!")
      }
    }

    if(!noBin) {
      temp <- binom2sumup_ages(parameters$Bin, parameters$Age)
      parameters <- c(parameters, Pos=list(temp$pos), Tot=list(temp$tot))
      parameters$Age <- temp$age
    }
  }

  res <- fun(parameters)

  message("Done.")
  res
}
