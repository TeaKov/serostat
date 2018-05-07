"pspline.fit"<-
function(response, x.var, ps.intervals = 8, wts = NULL, degree = 3, order = 3,
	link = "default", family = "gaussian", m.binomial = NULL, r.gamma =
	NULL, lambda = 0, x.predicted = NULL, ridge.adj = 0.0001)
{
# Function pspline.fit: univariate smoother using P-splines.
# Input: x.var= explanatory variable on abcissae.
# Input: response= response variable.
# Input: family=gaussian, binomial, poisson, Gamma distribution.
# Input: wts= vector of weights; default is vector of ones.
# Input: m.binomial=vector of binomial trials. Default is 1 vector.
# Input: r.gamma=vector of gamma shape parameters. Default is 1 vector.
# Input: link= link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical).
# Input: ps.intervals= number of intervals for B-splines. Default=8.
# Input: degree= degree of B-splines. Default=3.
# Input: order= order of difference penalty. Default=3.
# Input: lambda= smoothness regulalizing parameter ( >= 0). Default=0.
# Input: x.predicted=a list of x variables for prediction and twice stderr limits.
# Result: a scatterplot of (response, x.var) with smoothed fit and se bands.
# Output: A list: including AIC= deviance + 2*trace(Hat), dispers.parm, etc.
#
# Reference: Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and
#            penalties (with comments and rejoinder). Statistical Science, 11(2): 89-121.
#
#
# Support functions: pspline.checker(), pspline.fitter(), pspline.predictor()
#
#
# (c) 1995 Paul Eilers & Brian Marx
#
  y <- response
	x <- x.var
	if(base::missing(wts)) {
		wts <- base::rep(1, base::length(y))
	}
	parms <- pspline.checker(family, link, degree, order, ps.intervals,
		lambda, ridge.adj, wts)
	family <- parms$family
	link <- parms$link
	q <- parms$degree
	d <- parms$order
	ridge.adj <- parms$ridge.adj
	lambda <- parms$lambda
	ndx <- parms$ps.intervals
	wts <- parms$wts
	if(base::missing(m.binomial)) {
		m.binomial <- base::rep(1, base::length(y))
	}
	if(base::missing(r.gamma)) {
		r.gamma <- base::rep(1, base::length(y))
	}
	n <- base::length(y)
	xl <- base::min(x)
	xr <- base::max(x)
	xmax <- xr + 0.01 * (xr - xl)
	xmin <- xl - 0.01 * (xr - xl)
	dx <- (xmax - xmin)/ndx
	knots <- base::seq(xmin - q * dx, xmax + q * dx, by = dx)
	b <- splines::spline.des(knots, x, q + 1, 0 * x)$design
	n.col <- base::ncol(b)
	if(d < 0) {
		d <- base::min(3, (n.col - 1))
		warning(base::paste("penalty order cannot be negative: have used", d))
	}
	if((d - n.col + 1) > 0) {
		d <- n.col - 1
		warning(base::paste("penalty order was too large: have used", d))
	}
	if(ridge.adj > 0) {
		nix.ridge <- base::rep(0, n.col)
		p.ridge <- base::sqrt(ridge.adj) * base::diag(base::rep(1, n.col))
	}
	p <- base::diag(n.col)
	if(d != 0) {
		for(j in 1:d) {
			p <- base::diff(p)
		}
	}
	p <- base::sqrt(lambda) * p
	nix <- base::rep(0, n.col - d)
	b <- base::as.matrix(b)
	ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b,
		p, p.ridge, nix, nix.ridge, ridge.adj, wts)
	mu <- ps.fit$mu
	coef <- ps.fit$coef
	w <- ps.fit$w
	e <- 1e-009
	h <- stats::hat(ps.fit$f$qr, intercept = F)[1:n]
	trace <- base::sum(h) - 1
	if(family == "binomial") {
		dev <- 2 * base::sum((y + e) * base::log((y + e)/(mu + e)) + (m.binomial -
			y + e) * base::log((m.binomial - y + e)/(m.binomial - mu + e)
			))
		dispersion.parm <- 1
	}
	if(family == "poisson") {
		dev <- 2 * base::sum(y * log(y + e) - y - y * base::log(mu) + mu)
		dispersion.parm <- 1
	}
	if(family == "Gamma") {
		dev <- -2 * base::sum(r.gamma * (base::log((y + e)/mu) - ((y - mu)/mu)))
		ave.dev <- dev/n
		dispersion.parm <- (ave.dev * (6 + ave.dev))/(6 + 2 * ave.dev)
	}
	if(family == "gaussian") {
		dev <- base::sum(ps.fit$f$residuals^2)
		dispersion.parm <- dev/(n - trace)
	}
	aic <- dev + 2 * trace
	x.seq <- base::seq(xl, xr, length = 50)
	b.seq <- splines::spline.des(knots, x.seq, q + 1, 0 * x.seq)$design
	w.aug <- base::c(w, (nix + 1))
	yhat <- b.seq %*% base::as.vector(ps.fit$coef)
	half.meat <- base::sqrt(c(w)) * b
	meat <- base::t(half.meat) %*% half.meat
	if(ridge.adj > 0) {
		bread <- base::solve(meat + base::t(p) %*% p + base::t(p.ridge) %*% p.ridge)
	}
	if(ridge.adj == 0) {
		bread <- base::solve(meat + base::t(p) %*% p)
	}
	half.sw <- half.meat %*% bread
	var.beta <- base::t(half.sw) %*% half.sw
	var.yhat <- b.seq %*% var.beta %*% base::t(b.seq)
	stdev.yhat <- base::as.vector(base::sqrt(base::diag(var.yhat)))
	stdev.yhat <- base::sqrt(dispersion.parm) * stdev.yhat
	pivot <- 2 * stdev.yhat
	upper <- yhat + pivot
	lower <- yhat - pivot
	summary.yhat <- base::cbind(lower, yhat, upper)
	if(link == "logit") {
		summary.yhat <- 1/(1 + base::exp( - summary.yhat))
	}
	if(link == "probit") {
		summary.yhat <- base::apply(summary.yhat, base::c(1, 2), pnorm)
	}
	if(link == "cloglog") {
		summary.yhat <- (1 - base::exp( - base::exp(summary.yhat)))
	}
	if(link == "loglog") {
		summary.yhat <- base::exp( - base::exp( - summary.yhat))
	}
	if(link == "sqrt") {
		summary.yhat <- summary.yhat^2
	}
	if(link == "log") {
		summary.yhat <- base::exp(summary.yhat)
	}
	if(link == "recipical") {
		summary.yhat <- 1/(summary.yhat)
	}
	if(family == "binomial" && base::mean(m.binomial) != 1) {
		graphics::matplot(x.seq, summary.yhat, type = "l", lty = base::c(2, 1, 2), xlab
			 = "regressor", ylab = "estimated mean", main =
			"P-spline fit with twice std error bands")
	}
	if(mean(m.binomial) == 1) {
		graphics::matplot(x.seq, summary.yhat, type = "l", lty = c(2, 1, 2), xlab
			 = "regressor", ylab = "estimated mean", ylim = base::c(base::min(
			   base::min(y), base::min(summary.yhat[, 1])), base::max(max(y), base::max(
			summary.yhat[, 3]))), main =
			"P-spline fit with twice std error bands")
		graphics::matpoints(x, y, type = "p", pch = "O")
	}
	ps.predict <- NULL
	if(!base::missing(x.predicted)) {
		ps.predict <- pspline.predictor(x.predicted, knots, link, coef,
			q, var.beta, dispersion.parm)
	}
	llist <- base::list()
	llist$family <- family
	llist$link <- link
	llist$ps.intervals <- ndx
	llist$order <- d
	llist$degree <- q
	llist$lambda <- lambda
	llist$aic <- aic
	llist$deviance <- dev
	llist$eff.df <- trace
	llist$df.resid <- n - trace
	llist$dispersion.param <- dispersion.parm
	llist$summary.predicted <- ps.predict$summary.pred
	llist$coef <- coef
	llist
}

"pspline.fitter"<-
function(family, link, n.col, m.binomial, r.gamma, y, b, p, p.ridge, nix,
	nix.ridge, ridge.adj, wts)
{
	coef.est <- base::rep(1, base::ncol(b))
	if(family == "binomial") {
		mu <- (y + 0.5 * m.binomial)/2
	}
	if(family == "Gamma" || family == "poisson") {
		mu <- (y + 3)
	}
	if(family == "gaussian") {
		mu <- base::rep(base::mean(y), base::length(y))
	}
	it <- 0
	repeat {
		if(it == 0) {
			if(link == "identity") {
				eta <- mu
			}
			if(link == "log") {
				eta <- base::log(mu)
			}
			if(link == "sqrt") {
				eta <- base::sqrt(mu)
			}
			if(link == "logit") {
				eta <- base::log(mu/(m.binomial - mu))
			}
			if(link == "recipical") {
				eta <- 1/mu
			}
			if(link == "probit") {
				eta <- stats::qnorm(mu/m.binomial)
			}
			if(link == "cloglog") {
				eta <- base::log( - base::log(1 - mu/m.binomial))
			}
			if(link == "loglog") {
				eta <-  - base::log( - base::log(mu/m.binomial))
			}
		}
		it <- it + 1
		if(it > 25)
			break
		if(link == "identity") {
			mu <- eta
			h.prime <- 1
		}
		if(link == "log") {
			mu <- base::exp(eta)
			h.prime <- mu
		}
		if(link == "sqrt") {
			mu <- eta^2
			h.prime <- 2 * eta
		}
		if(link == "logit") {
			mu <- m.binomial/(1 + base::exp( - eta))
			h.prime <- mu * (1 - mu/m.binomial)
		}
		if(link == "recipical") {
			mu <- 1/eta
			h.prime <-  - (mu^2)
		}
		if(link == "probit") {
			mu <- m.binomial * stats::pnorm(eta)
			h.prime <- m.binomial * stats::dnorm(eta)
		}
		if(link == "cloglog") {
			mu <- m.binomial * (1 - base::exp( - base::exp(eta)))
			h.prime <- (m.binomial) * base::exp(eta) * base::exp( - base::exp(eta))
		}
		if(link == "loglog") {
			mu <- m.binomial * base::exp( - base::exp( - eta))
			h.prime <- m.binomial * base::exp( - eta) * base::exp( - base::exp( - eta
				))
		}
		if(family == "gaussian") {
			w <- base::rep(1, base::length(y))
		}
		if(family == "poisson") {
			w <- h.prime^2/mu
		}
		if(family == "binomial") {
			w <- h.prime^2/(mu * (1 - mu/m.binomial))
		}
		if(family == "Gamma") {
			w <- (r.gamma * h.prime^2)/mu^2
		}
		u <- (y - mu)/h.prime + eta
		if(ridge.adj > 0) {
			f <- lsfit(rbind(b, p, p.ridge), c(u, nix, nix.ridge),
				wt = c(wts, nix + 1, nix.ridge + 1) * c(w, (nix +
				1), (nix.ridge + 1)), intercept = F)
		}
		if(ridge.adj == 0) {
			f <- stats::lsfit(base::rbind(b, p), base::c(u, nix), wt = base::c(wts, nix + 1) *
			     base::c(w, (nix + 1)), intercept = F)
		}
		coef.old <- coef.est
		coef.est <- base::as.vector(f$coef)
		d.coef <- base::max(base::abs((coef.est - coef.old)/coef.old))
		if(d.coef < 1e-008)
			break
		print(base::c(it, d.coef))
		eta <- b %*% coef.est
	}
	if(it > 24) {
		warning(base::paste("parameter estimates did NOT converge in 25 iterations"
			))
	}
	llist <- base::list(coef = coef.est, mu = mu, f = f, w = w * wts)
	return(llist)
}

"pspline.checker"<-
function(family, link, degree, order, ps.intervals, lambda, ridge.adj, wts)
{
	if(link == "default" && family == "gaussian") {
		link <- "identity"
	}
	if(link == "default" && family == "poisson") {
		link <- "log"
	}
	if(link == "default" && family == "binomial") {
		link <- "logit"
	}
	if(link == "default" && family == "Gamma") {
		link <- "log"
	}
	if(family != "binomial" && family != "gaussian" && family != "poisson" &&
		family != "Gamma") {
		warning(base::paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"
			))
	}
	if((family == "binomial") && (link != "logit" && link != "probit" &&
		link != "cloglog" && link != "loglog")) {
		warning(base::paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"
			))
	}
	if((family == "Gamma") && (link != "log" && link != "recipical" && link !=
		"identity")) {
		warning(base::paste("Improper LINK option with family=Gamma. Choose: recipical, log, identity"
			))
	}
	if((family == "poisson") && (link != "log" && link != "sqrt" && link !=
		"identity")) {
		warning(base::paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"
			))
	}
	if((family == "gaussian") && (link != "identity")) {
		warning(base::paste("Improper LINK option with family=gaussian. Choose: identity"
			))
	}
	if(degree < 0) {
		degree <- 1
		warning(base::paste("degree must be non-neg integer: have used 1"))
	}
	if(order < 0) {
		order <- 0
		warning(base::paste("order must be non-neg integer: have used 0"))
	}
	if(ps.intervals < 2) {
		ps.intervals <- 2
		warning(base::paste("ps.intervals must be positive integer, > 1: have used 2"
			))
	}
	if(lambda < 0) {
		lambda <- 0
		warning(base::paste("lambda cannot be negative: have used 0"))
	}
	if(ridge.adj < 0) {
		ridge.adj <- 0
		warning(base::paste("ridge.adj cannot be negative: have used 0"))
	}
	if(min(wts) < 0) {
		warning(base::paste("At least one weight entry is negative"))
	}
	llist <- base::list(family = family, link = link, degree = degree, order =
		order, ps.intervals = ps.intervals, lambda = lambda, ridge.adj
		 = ridge.adj, wts = wts)
	return(llist)
}

"pspline.predictor"<-
function(x.predicted, knots, link, coef, q, var.beta, dispersion.parm)
{
	b.pred <- splines::spline.des(knots, x.predicted, q + 1, 0 * x.predicted)$design
	eta.pred <- b.pred %*% base::as.vector(coef)
	b.pred <- base::as.matrix(b.pred)
	if(base::length(x.predicted) > 1) {
		var.pred <- (b.pred) %*% var.beta %*% base::t(b.pred)
	}
	if(base::length(x.predicted) == 1) {
		var.pred <- base::t(b.pred) %*% var.beta %*% (b.pred)
	}
	stdev.pred <- base::as.vector(base::sqrt(base::diag(var.pred)))
	stdev.pred <- base::sqrt(dispersion.parm) * stdev.pred
	pivot <- base::as.vector(2 * stdev.pred)
	upper <- eta.pred + pivot
	lower <- eta.pred - pivot
	summary.pred <- base::cbind(lower, eta.pred, upper)
	if(link == "logit") {
		summary.pred <- 1/(1 + base::exp( - summary.pred))
	}
	if(link == "probit") {
		summary.pred <- base::apply(summary.pred, c(1, 2), pnorm)
	}
	if(link == "cloglog") {
		summary.pred <- (1 - base::exp( - base::exp(summary.pred)))
	}
	if(link == "loglog") {
		summary.pred <- base::exp( - base::exp( - summary.pred))
	}
	if(link == "sqrt") {
		summary.pred <- summary.pred^2
	}
	if(link == "log") {
		summary.pred <- base::exp(summary.pred)
	}
	if(link == "recipical") {
		summary.pred <- summary.predd <- 1/(summary.pred)
		summary.pred <- summary.predd[, 3:1]
	}
	summary.pred <- base::as.matrix(summary.pred)
	base::dimnames(summary.pred) <- base::list(NULL, base::c("-2std_Lower", "Predicted",
		"+2std_Upper"))
	llist <- base::list(summary.pred = summary.pred)
	return(llist)
}
