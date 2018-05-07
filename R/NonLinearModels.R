FarringtonABG <- function(d, p)
{
  result <- stats4::mle( start=base::list(alpha=p$alpha, beta=p$beta, gamma=p$gamma),
      method = "L-BFGS-B", lower = c(0.,0.,0.), upper = c(Inf, Inf, Inf),
      function(alpha, beta, gamma)
      {
        ab <- alpha/beta
        e_ba <- base::exp(-beta*d$Age)
        p <- 1-base::exp(ab*d$Age*e_ba+(1/beta)*(ab-gamma)*(e_ba-1)-gamma*d$Age)

        ll <- d$Pos*base::log(p)+(d$Tot-d$Pos)*base::log(1-p)
        return(-base::sum(ll))
      })
  result
}

FarringtonAB_G <- function(d, p)
{
  result <- stats4::mle( start=base::list(alpha=p$alpha, beta=p$beta), fixed=base::list(gamma=p$gamma),
      method = "L-BFGS-B", lower = c(0.,0.), upper = c(Inf, Inf),
      function(alpha, beta, gamma)
      {
        ab <- alpha/beta
        e_ba <- base::exp(-beta*d$Age)
        p <- 1-base::exp(ab*d$Age*e_ba+(1/beta)*(ab-gamma)*(e_ba-1)-gamma*d$Age)

        ll <- d$Pos*base::log(p)+(d$Tot-d$Pos)*base::log(1-p)
        return(-base::sum(ll))
      })
  result
}

NonLinearModels <- function(data, parameters, fun)
{
  UseMethod("NonLinearModels", fun)
}

NonLinearModels.default <- function(data, parameters, fun)
{
  noTot <- !("Tot" %in% base::names(data))
  noPos <- !("Pos" %in% base::names(data))
  noAge <- !("Age" %in% base::names(data))

  noAlpha <- !("alpha" %in% base::names(parameters))
  noBeta  <- !("beta"  %in% base::names(parameters))
  noGamma <- !("gamma" %in% base::names(parameters))

  if(noTot)
    warning("Tot not specified!")
  if(noPos)
    warning("Pos not specified!")
  if(noAge)
    warning("Age not specified!")
  if(noAlpha)
    warning("Alpha not specified!")
  if(noBeta)
    warning("Beta not specified!")
  if(noGamma)
    warning("Gamma not specified!")
  if(noTot || noPos || noAge)
    stop("Please specify the necessary parameters!")

  res <- fun(data, parameters)

  # class(res) <- c("NonLinearModels", as.character(substitute(fun)))
  message("Done. Pass result to summary() for MLE results.")
  res
}

# print.NonLinearModels <- function(x, ...)
# {
#   print(paste("Model: ", class(x)[2]), quote = F)
#   classy <- "mle"; attr(classy, "package") <- "stats4";
#   class(x) <- classy
#   summary(x)
# }
