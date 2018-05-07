RAS_FOI <- function(t, states, params)
{
  n.ages <- params[1]
  mu <- params[2:(n.ages+1)]
  foi <- params[(n.ages+2):(2*n.ages+1)]
  alpha.rate <- params[(2*n.ages+2)]

  S <- states[           1:(n.ages)]
  I <- states[  (n.ages+1):(2*n.ages)]
  R <- states[(2*n.ages+1):(3*n.ages)]

  dS <- - (foi+mu)*S                # Susceptibles
  dI <- +foi*S -(alpha.rate+mu)*I   # Infection
  dR <- +alpha.rate*I-mu*R          # Immune

  list(c(dS,dI,dR))
}

RAS_BETAS = function(t, states, params)
{
  n.ages <- params[1]
  mu <- params[2:(n.ages+1)]
  betas <- matrix(params[(n.ages+2):((n.ages+1)*n.ages+1)],n.ages,n.ages)
  alpha.rate <- params[((n.ages+1)*n.ages+2)]

  S <- states[           1:(n.ages)]
  I <- states[  (n.ages+1):(2*n.ages)]
  R <- states[(2*n.ages+1):(3*n.ages)]

  lambda <- apply(betas*I,2,"sum")

  dS <- -(lambda+mu)*S                # Susceptibles
  dI <- +lambda*S -(alpha.rate+mu)*I  # Infection
  dR <- +alpha.rate*I-mu*R            # Immune

  list(c(dS,dI,dR))
}

RAS_PDE = function(t, states, params)
{
  alpha <- 1; # dt/dage
  n.ages <- params[1]
  mu <- params[2:(n.ages+1)]
  betas <- matrix(params[(n.ages+2):((n.ages+1)*n.ages+1)],n.ages,n.ages)
  alpha.rate <- params[((n.ages+1)*n.ages+2)]

  S <- states[           1:(n.ages)]
  I <- states[  (n.ages+1):(2*n.ages)]
  R <- states[(2*n.ages+1):(3*n.ages)]

  lambda = apply(betas*I,2,"sum")

  dS <- -(lambda+mu)*S - alpha * S + alpha * c(0, S[1:(n.ages-1)])                # Susceptible
  dI <- +lambda*S -(alpha.rate+mu)*I - alpha * I + alpha * c(0, I[1:(n.ages-1)])  # Infected
  dR <- +alpha.rate*I-mu*R - alpha * R + alpha * c(0, R[1:(n.ages-1)])            # Recovered

  dS[1] <- dS[1] - (sum(dS) + sum(dI) + sum(dR))

  list(c(dS,dI,dR))
}

solveRAS <- function(config, parameters, state, fun, ...)
{
  UseMethod("solveRAS", fun)
}

solveRAS.default <- function(config, parameters, state, fun, ...)
{
  function_string = base::as.character(substitute(fun))
  noPDE <- FALSE

  noTinit <- !("Tinit" %in% base::names(config))
  noTruns <- !("Truns" %in% base::names(config))
  noTol <- !("tol" %in% base::names(config))
  noRes <- !("res" %in% base::names(config))

  if(noTinit)
    warning("Tinit not specified!")
  if(noTruns)
    warning("Truns not specified!")
  if(noTol)
    warning("tol not specified!")
  if(noRes)
    warning("res not specified!")
  if(noTinit || noTruns || noTol || noRes)
    stop("Please specify the necessary parameters!")

  if(function_string == "RAS_FOI" || function_string == "RAS_BETAS") {
    noPDE <- TRUE
    noCohortSize <- !("cohort.size" %in% base::names(config))
    if(noCohortSize) {
      warning("cohort.size not specified!")
      stop("Please specify the necessary parameters!")
    }
  }

  fullout <- base::array(NA,dim=base::c(config$Tinit[2],config$res,300))

  for (initrun in config$Tinit[1]:config$Tinit[2])
  {
    # Some starting parameters
    print(paste("Initialisation run : Year = ",as.character(initrun),sep=""))
    states = ifelse(states < config$tol, 0, states)
    # moving everyone, one state forward
    if (noPDE && initrun != config$Tinit[1]) {
      for (j in 99:1)
      {
        states[j+  1] <- states[j    ]
        states[j+101] <- states[j+100]
        states[j+201] <- states[j+200]
      }
      states[1]   <- config$cohort.size 	# Completely susceptible at birth
      states[101] <- 0
      states[201] <- 0
    }

    # Time steps determined by resolution
    times  <- base::seq(0, 1, length=config$res)

    # Output of the system of ODEs
    out <- base::as.data.frame(deSolve::lsoda(states, times, fun, parameters, hmax=1.0, ...)[,-1])
    states <- base::as.matrix(out[config$res,])

    fullout[initrun,,] <- base::as.matrix(out)
  }

  res <- list(fullout=fullout, out=out)
}
