SmoothSplines <- function(p)
{
  gam::gam(p$y~gam::s(p$a, df=p$dfopts), family=binomial(link=p$link))
}

CubicRegSplines <- function(p)
{
  mgcv::gam(p$y~s(p$a,bs="cr"),family=binomial(link=p$link))
}

ThinPlateRegSplines <- function(p)
{
  mgcv::gam(p$y~s(p$a,bs="tp"),family=binomial(link=p$link))
}

CubicRegSplinesMM <- function(p)
{
  mgcv::gamm(p$y~s(p$a,bs="cr"),family=binomial(link=p$link))
}

ThinPlateRegSplinesMM <- function(p)
{
  mgcv::gamm(p$y~s(p$a,bs="tp"),family=binomial(link=p$link))
}

AdaptiveSplineSmoothing <- function(p)
{
  with(p, {
    knots_l <<- knots
    if(adap)
    {
      var.knot_l <<- var.knot
      AdaptFit::asp(y~f(a,knots=knots_l,var.knot=var.knot_l), tol=1e-3, niter=1000, niter.var=100000, omit.missing=TRUE, adap=TRUE, family="binomial", spar.method="REML")
    }
    else
      AdaptFit::asp(y~f(a, knots=knots_l), adap=FALSE, family="binomial", tol=1e-6, niter=1000, niter.var=100000, omit.missing=TRUE, spar.method="REML")
  })
}

SemiParametricRegression <- function(parameters, fun)
{
  UseMethod("SemiParametricRegression", fun)
}

SemiParametricRegression.default <- function(parameters, fun)
{
  fs = as.character(substitute(fun))

  noY <- !("y" %in% names(parameters))
  noA <- !("a" %in% names(parameters))

  unspecifiedParams <- F

  if(fs == "SmoothSplines")
  {
    if(!("dfopts" %in% names(parameters))) {
      warning("dfopts not specified!")
      unspecifiedParams <- T
    }

    if(!("link" %in% names(parameters))) {
      warning("link not specified!")
      unspecifiedParams <- T
    }
  }

  if(fs == "CubicRegSplines" || fs == "ThinPlateRegSplines" || fs == "CubicRegSplinesMM" || fs == "ThinPlateRegSplinesMM")
  {
    if(!("link" %in% names(parameters))) {
      warning("link not specified!")
      unspecifiedParams <- T
    }
  }

  if(fs == "AdaptiveSplineSmoothing")
  {
    if(!("adap" %in% names(parameters)))
    {
      warning("adap not specified!")
      unspecifiedParams <- T
    }
    else
    {
      if(!("knots" %in% names(parameters)))
      {
        warning("knots not specified!")
        unspecifiedParams <- T
      }
      if(parameters$adap && !("var.knot" %in% names(parameters)))
      {
        warning("var.knot not specified!")
        unspecifiedParams <- T
      }
    }
  }

  if(noY) {
    warning("y not specified!")
    unspecifiedParams <- T
  }
  if(noA) {
    warning("a not specified!")
    unspecifiedParams <- T
  }
  if(unspecifiedParams)
    stop("Please specify the necessary parameters!")

  res <- fun(parameters)
  message("Done.")
  res
}

