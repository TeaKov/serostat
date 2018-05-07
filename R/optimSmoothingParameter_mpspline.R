optimSmoothingParameter <- function(mpspline.params, alpha.init=2, useAIC=TRUE, upper=100) {
  if(useAIC)
    message("Using AIC for optimization")
  else
    message("Using BIC for optimization")

  noResponse <- !("response" %in% base::names(mpspline.params))
  noXVar <- !("x.var" %in% base::names(mpspline.params))

  if(noResponse)
    warning("response not specified!")
  if(noXVar)
    warning("x.var not specified!")
  if(noResponse || noXVar)
    stop("Please specify the necessary parameters!")

  m.response <- mpspline.params$response
  m.x.var <- mpspline.params$x.var

  m.ps.intervals <- 20
  if("ps.intervals" %in% base::names(mpspline.params))
    m.ps.intervals <- mpspline.params$ps.intervals

  m.degree <- 3
  if("degree" %in% base::names(mpspline.params))
    m.degree <- mpspline.params$degree

  m.order <- 2
  if("order" %in% base::names(mpspline.params))
    m.order <- mpspline.params$order

  m.link <- "identity"
  if("link" %in% base::names(mpspline.params))
    m.link <- mpspline.params$link

  m.family <- "gaussian"
  if("family" %in% base::names(mpspline.params))
    m.family <- mpspline.params$family

  m.kappa <- 1e8
  if("kappa" %in% base::names(mpspline.params))
    m.kappa <- mpspline.params$kappa

  res <- optim(par = alpha.init, method="Brent", fn = function(x) {
    fit <- mpspline.fit(response=m.response, x.var=m.x.var, ps.intervals=m.ps.intervals,
                        degree=m.degree, order=m.order, link=m.link, family=m.family, alpha=x, kappa=m.kappa)
    if(useAIC)
      res <- fit$aic
    else
      res <- fit$bic
    res
  }, lower=0, upper=upper)

  res
}
