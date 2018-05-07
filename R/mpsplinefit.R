# Monotone logistic B-spline fit with logit-, cloglog- and probit-link function

# Function
mpspline.fit<-function(response, x.var, ps.intervals=20, degree=3, order=2, link="identity", family="gaussian", alpha=2, kappa=1e8)
{
	y <- response
	x <- x.var
	wts <- base::rep(1, base::length(y))
	q <- degree
	d <- order
	ndx <- ps.intervals
	m.binomial <- base::rep(1, base::length(y))
	n <- base::length(y)
	xl <- base::min(x)
	xr <- base::max(x)
	xmax <- xr + 0.01 * (xr - xl)
	xmin <- xl - 0.01 * (xr - xl)
	dx <- (xmax - xmin)/ndx
	knots <- base::seq(xmin - q * dx, xmax + q * dx, by = dx)
	b <- splines::spline.des(knots, x, q + 1, 0 * x)$design
	n.col <- base::ncol(b)
	p <- base::sqrt(alpha)*base::diff(base::diff(base::diag(n.col)))
	mp <- base::sqrt(kappa)*base::diff(base::diag(n.col))
	nix <- base::rep(0, n.col - d)
	mnix <- base::rep(0, n.col - 1)
	b <- base::as.matrix(b)
	coef.est <- base::rep(1, ncol(b))
	ineqdiag <- base::diag(base::diff(coef.est)>=1e-9)

	if(family == "gaussian") {
		mu <- base::rep(mean(y),base::length(y))
	}

	if(family == "binomial") {
		mu <- (y + 0.5 * m.binomial)/2
	}

	it <- 0
	repeat {
		if(it == 0) {
			if(link == "identity") {
				eta <- mu
			}
			if(link == "logit") {
				eta <- base::log(mu/(m.binomial - mu))
			}
			if(link == "probit") {
				eta <- stats::qnorm(mu/m.binomial)
			}
			if(link == "cloglog") {
				eta <- base::log( - base::log(1 - mu/m.binomial))
			}
		}

		it <- it + 1
		if(it > 25)
			break
		if(link == "identity") {
			mu <- eta
			h.prime <- 1
		}

		if(link == "logit") {
			mu <- m.binomial/(1 + base::exp( - eta))
			h.prime <- mu * (1 - mu/m.binomial)
		}

		if(link == "probit") {
			mu <- m.binomial * stats::pnorm(eta)
			h.prime <- m.binomial * stats::dnorm(eta)
		}

		if(link == "cloglog") {
			mu <- m.binomial * (1 - base::exp( - base::exp(eta)))
			h.prime <- (m.binomial) * base::exp(eta) * base::exp( - base::exp(eta))
		}

		if(family == "gaussian") {
			w <- base::rep(1,base::length(y))
		}

		if(family == "binomial") {
			w <- h.prime^2/(mu * (1 - mu/m.binomial))
		}

		u <- (y - mu)/h.prime + eta
		f <- stats::lsfit(base::rbind(b, p, ineqdiag%*%mp), base::c(u, nix, mnix), wt = base::c(wts, nix + 1, mnix + 1) *
		  base::c(w, (nix + 1), (mnix + 1)), intercept = F)
		coef.old <- coef.est
		coef.est <- base::as.vector(f$coef)
		ineqdiag<-base::diag(base::diff(coef.est)<0)
		d.coef <- base::max(base::abs((coef.est - coef.old)/coef.old))

		if(d.coef < 1e-20)
			break
		#print(c(it, d.coef))
		eta <- b %*% coef.old
	}

	w <- w*wts
	e <- 1e-009
	h <- stats::hat(f$qr, intercept = F)[1:n]
	trace <- base::sum(h) - 1

	if(family == "gaussian") {
		dev <- base::sum((y-eta)^2)
		dispersion.parm <- dev/(n - trace)
	}

	if(family == "binomial") {
		dev <- 2 * base::sum((y + e) * base::log((y + e)/(mu + e)) + (m.binomial -
			y + e) * base::log((m.binomial - y + e)/(m.binomial - mu + e)))
		dispersion.parm <- 1
	}

	aic <- dev + 2 * trace
	bic <- dev + base::log(n) * trace
	x.seq <- base::seq(xl, xr, by=1)
	b.seq <- splines::spline.des(knots, x.seq, q + 1, 0 * x.seq)$design
	w.aug <- base::c(w, (nix + 1))
	yhat <- b.seq %*% base::as.vector(coef.old)
	summary.yhat <- yhat

	if(link == "logit") {
		summary.yhat <- 1/(1 + base::exp( - summary.yhat))
	}

	if(link == "probit") {
		summary.yhat <- base::apply(summary.yhat, c(1, 2), stats::pnorm)
	}

	if(link == "cloglog") {
		summary.yhat <- (1 - base::exp( - base::exp(summary.yhat)))
	}

	return(base::list(x=x.seq,yhat=summary.yhat,aic=aic,bic=bic,dev=dev))
}
