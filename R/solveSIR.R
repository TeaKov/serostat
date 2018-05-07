# This contains the simple SIR examples from Chapter 3

SIR_NoVac <- function(t, state, parameters)
{
  with(base::as.list(base::c(state,parameters)),
  {
    dS = N*mu - beta*I*S - mu*S
    dI = beta*I*S - (nu+alpha+mu)*I
    dR = nu*I - mu*R
    base::list(c(dS, dI, dR))
  })
}

SIR_Vac<- function(t, state, parameters)
{
  with(base::as.list(base::c(state,parameters)),
  {
    dS = N*mu*(1-p) - beta*I*S - mu*S
    dI = beta*I*S - (nu+alpha+mu)*I
    dR = nu*I - mu*R + N*mu*p
    base::list(c(dS, dI, dR))
  })
}

# Normalized values!
# s = S/N
# i = I/N
# r = R/N
SIR_ConstFOI<-function(t,state,parameters)
{
  with(base::as.list(base::c(state, parameters)),
  {
    ds <- -lambda*s
    di <- lambda*s - nu*i
    dr <- nu*r
    base::list(c(ds, di, dr))
  })
}

SIR_SubPop <- function(t,state,parameters)
{
  with(base::as.list(base::c(state, parameters)),
  {
    betatilde_m <- base::matrix(params[1:(ndims*ndims)],ncol=ndims, nrow=ndims)
    nu_m <- params[((ndims*ndims)+1):((ndims*ndims)+nparams)]

    s_m <- state[1:nparams]
    i_m <- state[(nparams+1):(2*nparams)]
    r_m <- state[((2*nparams)+1):(3*nparams)]

    temp <- t(t(betatilde_m) %*% i_m)
    ds <- -(temp) * s_m + mu - mu * s_m
    di <- temp * s_m - nu_m * i_m - mu * i_m
    dr <- nu_m * i_m - mu * r_m
    base::list(base::c(ds,di,dr))
  })
}

solveSIR <- function(t, parameters, state, fun)
{
  UseMethod("solveSIR", fun)
}

solveSIR.default <- function(t, parameters, state, fun)
{
  function_string = base::as.character(substitute(fun))

  noNumParams <- !("nparams" %in% base::names(parameters))
  noNumDims <- !("ndims" %in% base::names(parameters))
  if(function_string == "SIR_SubPop" && noNumParams)
    warning("nparams not specified!")
  if(function_string == "SIR_SubPop" && noNumDims)
    warning("ndims not specified!")
  if(function_string == "SIR_SubPop" && (noNumParams || noNumDims))
    stop("Please specify necessary parameters!")

  res <- base::as.data.frame(deSolve::ode(y=state, times=t, func=fun, parms=parameters))
  res
}
