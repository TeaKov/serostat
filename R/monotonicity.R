pavit <- function(pos=pos,tot=rep::rep(1,base::length(pos)))
{
  gi <- pos/tot
  pai1 <- pai2 <- gi
  N <- base::length(pai1)
  ni <- tot
  for(i in 1:(N - 1)) {
    if(pai2[i] > pai2[i + 1]) {
      pool <- (ni[i]*pai1[i] + ni[i+1]*pai1[i + 1])/(ni[i]+ni[i+1])
      pai2[i:(i + 1)] <- pool
      k <- i + 1
      for(j in (k - 1):1) {
        if(pai2[j] > pai2[k]) {
          pool.2 <- base::sum(ni[j:k]*pai1[j:k])/(base::sum(ni[j:k]))
          pai2[j:k] <- pool.2
        }
      }
    }
  }
  return(base::list(pai1=pai1,pai2=pai2))
}

foi.num<-function(x,p)
{
  grid <- base::sort(base::unique(x))
  pgrid <- (p[base::order(x)])[base::duplicated(base::sort(x))==FALSE]
  dp <- base::diff(pgrid)/base::diff(grid)
  foi <- stats::approx((grid[-1]+grid[-base::length(grid)])/2,dp,grid[c(-1,-base::length(grid))])$y/(1-pgrid[c(-1,-base::length(grid))])
  return(base::list(grid=grid[c(-1,-base::length(grid))],foi=foi))
}

pwcfoi <- function(response, x.var, breaks){
  h1 <- function(beta1)
  {
    n <- base::length(response)
    int <- base::rep(0,base::length(x.var))
    integrand <- base::rep(0,base::length(x.var))

    for (i in 1:(base::length(breaks)-1)) {
      int <- int+((breaks[i+1]-breaks[i])*beta1[i]*(x.var>=breaks[i+1])+(x.var-breaks[i])*beta1[i]*(x.var<breaks[i+1])*(x.var>breaks[i]))
      integrand <- integrand+beta1[i]*(x.var<breaks[i+1])*(x.var>=breaks[i])
    }

    lambda <- integrand
    clambda <- int
    pihat <- 1-base::exp(-clambda)
    pli <- rep(0,n)
    pli[response==0] <- (1-pihat)[response==0]
    pli[response==1] <- pihat[response==1]
    lpli <- base::log(pli+1e-6)

    return(base::list(pli=-base::sum(lpli),x.var=x.var,prob=pihat,foi=lambda))
  }

  h2 <- function(beta) { return(h1(beta)$pli) }

  result.nlm <- stats::nlm(f=h2,p=base::rep(0.1,base::length(breaks)-1),iterlim=500,print.level=0)
  prev <- h1(result.nlm$estimate)$prob
  foi <- h1(result.nlm$estimate)$foi
  dev <- 2*h1(result.nlm$estimate)$pli
  k <- base::length(breaks)-1
  aic <- dev+2*k
  bic <- dev+base::log(base::length(response))*k

  return(base::list(prev=prev,foi=foi,lambda.vec=result.nlm$estimate,x.var=x.var,aic=aic,bic=bic,dev=dev))
}

pwcrate <- function(y.var, x.var, n.var=rep(1,length(y.var)), breaks, startpar=rep(2e-1,length(breaks)-1)) {
  h1 <- function(parms) {
    parms <- parms^2
    n <- base::length(y.var)
    int <- rep(0, base::length(x.var))
    integrand <- base::rep(0,length(x.var))

    for (i in 1:(base::length(breaks)-1)) {
      int <- int+((breaks[i+1]-breaks[i])*parms[i]*(x.var>=breaks[i+1])+(x.var-breaks[i])*parms[i]*(x.var<breaks[i+1])*(x.var>breaks[i]))
      integrand <- integrand+parms[i]*(x.var<breaks[i+1])*(x.var>=breaks[i])
    }

    rate <- integrand
    cumrate <- int
    probhat <- 1-base::exp(-cumrate)
    lli <- y.var*base::log(probhat+1e-12)+(n.var-y.var)*base::log(1-probhat+1e-12)

    return(base::list(ll=-2*base::sum(lli),x.var=x.var,prob=probhat,rate=rate))
  }

  h2 <- function(parms) { return(h1(parms)$ll) }

  result.nlm <- stats::nlm(f=h2,p=startpar,iterlim=500,print.level=0,hessian=TRUE)
  se <- base::sqrt(base::diag(solve(result.nlm$hessian)))
  result <- h1(result.nlm$estimate)
  k <- base::length(breaks)-1
  dev <- result$ll
  aic <- dev+2*k
  bic <- dev+base::log(base::sum(n.var))*k

  return(base::list(prob=result$prob,rate=result$rate^2,aic=aic,bic=bic,deviance=dev,ratevec=base::round(result.nlm$estimate^2,8),breaks=breaks,x.var,y.var,n.var))
}

KeidingMonotonicity <- function(p)
{
  xx <- pavit(pos=p$pos, tot=p$tot)
  foi.k <- foi.num(p$grid, xx$pai2)$foi
  foi.k[base::is.na(foi.k)] <- 0
  foi.k[foi.k>10] <- 0
  age.k <- foi.num(p$grid, xx$pai2)$grid
  fit.k<- stats::ksmooth(age.k ,foi.k ,kernel="normal", bandwidth=p$bw, n.points=base::length(age.k))
  pihat <- 1-base::exp(-base::cumsum(c(age.k[1], base::diff(age.k))))*fit.k$y
}

