GF <- function(PP, PN, NP, NN, a, alpha1eta=0.06, beta1eta=0.2, gamma1eta=0.005,
               alpha2eta=0.06, beta2eta=0.2, gamma2eta=0.005, thetaeta=1) {
  alpha1 <- base::exp(alpha1eta);
  beta1 <- base::exp(beta1eta);
  gamma1 <- base::exp(gamma1eta);
  alpha2 <- base::exp(alpha2eta);
  beta2 <- base::exp(beta2eta);
  gamma2 <- base::exp(gamma2eta);
  theta <- base::exp(thetaeta);
  Lambda1 <- alpha1*gamma1^(beta1+1)*base::gamma(beta1+1)*stats::pgamma(a/gamma1,beta1+1)
  Lambda2 <- alpha2*gamma2^(beta2+1)*base::gamma(beta2+1)*stats::pgamma(a/gamma2,beta2+1)
  p00 <- (base::exp(Lambda1/theta)+base::exp(Lambda2/theta)-1)^(-theta)
  p10 <- base::exp(-Lambda2)-p00
  p01 <- base::exp(-Lambda1)-p00
  p11 <- 1-p00-p01-p10
  return(-base::sum(PP*base::log(p11)+PN*base::log(p10)+NP*base::log(p01)+NN*base::log(p00)))
}
