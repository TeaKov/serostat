estimateLifeExpectancy <- function(nd, ps, age) {
  demfit<-mgcv::gam(nd~s(age),offset=base::log(ps),family="poisson",link="log")
  muy<-stats::predict(demfit,type="response")
  My<-base::exp(-base::cumsum(muy))
  L<-base::mean(My)*100
  list(LifeExpectancy=L, My=My, muy=muy)
}


waifw.6parms <- function(foihat, muhat, breaks, N, D, Lmax, muy) {
  L <- Lmax * base::mean(base::exp(-base::cumsum(muhat*(base::floor(breaks)[-1]-base::floor(breaks)[-base::length(breaks)]))))
  phi <- base::cumsum((foihat+muhat)*(breaks[-1]-breaks[-7]))
  psi <- N*D/L*foihat/(foihat+muhat)*(base::exp(-c(0,phi[-6]))-base::exp(-phi))
  psi[foihat+muhat==0]<-N*D/L*(base::exp(-c(0,phi[-6]))-base::exp(-phi))[foihat+muhat==0]
  My <- base::exp(-base::cumsum(muy))

  # W1
  ######
  msg1 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6)
  R0vec <- NA
  bij <- NA
  Dij <- t(matrix(c(
    psi[1],0,0,0,0,sum(psi[-1]),
    0,psi[2],0,0,0,sum(psi[-2]),
    0,0,psi[3],0,0,sum(psi[-3]),
    0,0,0,psi[4],0,sum(psi[-4]),
    0,0,0,0,psi[5],sum(psi[-5]),
    0,0,0,0,0,sum(psi[1:6])),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg1<-c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij <- matrix(beta[6],nrow=Lmax,ncol=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- beta[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- beta[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- beta[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- beta[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- beta[6]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
    }
  }

  if (det(Dij)==0)
    msg1 <- base::c("D(lambda) is not invertible")
  beta1 <- beta
  R0vec1 <- R0vec
  bij1 <- bij

  # W2
  ######
  msg2 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6)
  R0vec <- NA
  bij <- NA
  Dij <- base::t(matrix(c(
    base::sum(psi[1:2]),0,psi[3],psi[4],psi[5],psi[6],
    psi[1],psi[2],psi[3],psi[4],psi[5],psi[6],
    0,0,base::sum(psi[1:3]),psi[4],psi[5],psi[6],
    0,0,0,base::sum(psi[1:4]),psi[5],psi[6],
    0,0,0,0,base::sum(psi[1:5]),psi[6],
    0,0,0,0,0,base::sum(psi[1:6])),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg2 <- base::c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij<-matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[3],1:breaks[2]] <- beta[1]
      bij[1:breaks[2],(breaks[2]+1):breaks[3]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- beta[2]
      bij[1:breaks[4],(breaks[3]+1):breaks[4]] <- beta[3]
      bij[(breaks[3]+1):breaks[4],1:breaks[4]] <- beta[3]
      bij[1:breaks[5],(breaks[4]+1):breaks[5]] <- beta[4]
      bij[(breaks[4]+1):breaks[5],1:breaks[5]] <- beta[4]
      bij[1:breaks[6],(breaks[5]+1):breaks[6]] <- beta[5]
      bij[(breaks[5]+1):breaks[6],1:breaks[6]] <- beta[5]
      bij[1:breaks[7],(breaks[6]+1):breaks[7]] <- beta[6]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- beta[6]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
    }
  }
  if (det(Dij)==0)
    msg2 <- base::c("D(lambda) is not invertible")
  beta2 <- beta
  R0vec2 <- R0vec
  bij2 <- bij

  # W3
  ######
  msg3 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6);R0vec<-NA;bij<-NA
  Dij <- base::t(matrix(base::c(
    base::sum(psi[1:3]),0,0,psi[4],psi[5],psi[6],
    psi[1],psi[2],psi[3],psi[4],psi[5],psi[6],
    psi[1],0,base::sum(psi[2:3]),psi[4],psi[5],psi[6],
    0,0,0,base::sum(psi[1:4]),psi[5],psi[6],
    0,0,0,0,base::sum(psi[1:5]),psi[6],
    0,0,0,0,0,base::sum(psi[1:6])),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg3<- base::c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij <- matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[4],1:breaks[2]] <- beta[1]
      bij[1:breaks[2],1:breaks[4]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- beta[2]
      bij[(breaks[2]+1):breaks[4],(breaks[3]+1):breaks[4]] <- beta[3]
      bij[(breaks[3]+1):breaks[4],(breaks[2]+1):breaks[4]] <- beta[3]
      bij[1:breaks[5],(breaks[4]+1):breaks[5]] <- beta[4]
      bij[(breaks[4]+1):breaks[5],1:breaks[5]] <- beta[4]
      bij[1:breaks[6],(breaks[5]+1):breaks[6]] <- beta[5]
      bij[(breaks[5]+1):breaks[6],1:breaks[6]] <- beta[5]
      bij[1:breaks[7],(breaks[6]+1):breaks[7]] <- beta[6]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- beta[6]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(base::c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
    }
  }

  if (det(Dij)==0)
    msg3 <- base::c("D(lambda) is not invertible")
  beta3 <- beta
  R0vec3 <- R0vec
  bij3 <- bij

  # W4
  ######
  msg4 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6)
  R0vec <- NA
  bij <- NA
  Dij <- base::t(matrix(base::c(
    base::sum(psi[1:6]),0,0,0,0,0,
    0,base::sum(psi[1:6]),0,0,0,0,
    0,0,base::sum(psi[1:6]),0,0,0,
    0,0,0,base::sum(psi[1:6]),0,0,
    0,0,0,0,base::sum(psi[1:6]),0,
    0,0,0,0,0,base::sum(psi[1:6])),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg4 <- base::c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij<-matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[7]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],1:breaks[7]] <- beta[2]
      bij[(breaks[3]+1):breaks[4],1:breaks[7]] <- beta[3]
      bij[(breaks[4]+1):breaks[5],1:breaks[7]] <- beta[4]
      bij[(breaks[5]+1):breaks[6],1:breaks[7]] <- beta[5]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- beta[6]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(base::c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
    }
  }
  if (det(Dij)==0)
    msg4 <- base::c("D(lambda) is not invertible")
  beta4 <- beta
  R0vec4 <- R0vec
  bij4 <- bij

  # W5
  ######
  msg5 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6)
  R0vec <- NA
  bij <- NA
  Dij <- base::t(matrix(base::c(
    psi[1],0,0,0,0,base::sum(psi[-1]),
    0,psi[2],0,0,0,base::sum(psi[-2]),
    0,0,psi[3],0,0,base::sum(psi[-3]),
    0,0,0,psi[4],0,base::sum(psi[-4]),
    0,0,0,0,psi[5],base::sum(psi[-5]),
    0,0,0,0,psi[6],base::sum(psi[1:5])),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg5 <- base::c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij <- matrix(beta[6],ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- beta[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- beta[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- beta[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- beta[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- beta[5]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(base::c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
    }
  }
  if (det(Dij)==0)
    msg5 <- base::c("D(lambda) is not invertible")
  beta5 <- beta
  R0vec5 <- R0vec
  bij5 <- bij

  # W6
  ######
  msg6 <- base::c("Regular Configuration")
  beta <- base::rep(NA,6)
  R0vec <- NA
  bij <- NA
  Dij <- base::t(matrix(base::c(
    psi[1],0,0,0,0,0,
    0,psi[2],0,0,0,0,
    0,0,psi[3],0,0,0,
    0,0,0,psi[4],0,0,
    0,0,0,0,psi[5],0,
    0,0,0,0,0,psi[6]),ncol=6))

  if(det(Dij)!=0) {
    beta <- solve(Dij,foihat)
    if(sum(beta<0)>0)
      msg6 <- base::c("D(lambda)lambda is irregular")
    if(sum(beta<0)==0) {
      bij <- matrix(0,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- beta[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- beta[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- beta[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- beta[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- beta[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- beta[6]
      R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
      Mij <- base::diag(base::c(My))
      R0vec <- base::as.double(base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values)
      print(base::c(base::round(base::max(R0vec),3)))
    }
  }
  if (det(Dij)==0)
    msg6 <- base::c("D(lambda) is not invertible")
  beta6 <- beta
  R0vec6 <- R0vec
  bij6 <- bij
  return(base::list(w1=base::list(message=msg1,beta=beta1,R0hat=base::max(R0vec1),bij=bij1,L=L,N=N,D=D),w2=base::list(message=msg2,beta=beta2,R0hat=base::max(R0vec2),bij=bij2,L=L,N=N,D=D),w3=base::list(message=msg3,beta=beta3,R0hat=base::max(R0vec3),bij=bij3,L=L,N=N,D=D),
                    w4=base::list(message=msg4,beta=beta4,R0hat=base::max(R0vec4),bij=bij4,L=L,N=N,D=D),w5=base::list(message=msg5,beta=beta5,R0hat=base::max(R0vec5),bij=bij5,L=L,N=N,D=D),w6=base::list(message=msg6,beta=beta6,R0hat=base::max(R0vec6),bij=bij6,L=L,N=N,D=D)))
}

waifw.fitter <- function(a, y, n=rep(1,length(y)), waifw.choice, muy, breaks, N, D, Lmax, A, startpar=NULL, plots="TRUE") {
  L <- Lmax* base::mean(base::exp(-base::cumsum(muy)))
  My <- base::exp(-base::cumsum(muy))
  
  # Monotonized piecewise constant FOI
  pcw.fit <- pwcrate(y.var=y, x.var=a, breaks=breaks)
  foihat <- pcw.fit$ratevec
  
  # Constant mortality rate
  muhat <- rep(NA,(length(breaks)-1))
  for (i in 1:(length(breaks)-1))
    muhat[i] <- mean(muy[(floor(breaks)[i]+1):(floor(breaks)[i+1])])
  
  # The result of the following function will be that W1, and W6 are invertible
  # W2, and W3 are irregular, and W4, and W5 are the only regulars
  waifw6.fit <- waifw.6parms(foihat=foihat, muhat=muhat, breaks=breaks, N=N, D=D, Lmax=Lmax, muy=muy)

  waifwproc <- function(a,y,b,n=rep(1,length(y)),waifw=waifw.choice,muy,breaks,Lmax,N,D,plots="TRUE") {
    b <- base::round(b^2,8)
    if (waifw==1) {
      bij <- matrix(b[6],nrow=Lmax,ncol=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- b[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- b[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- b[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- b[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- b[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- b[6]
    }

    if (waifw==2) {
      bij <- matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[3],1:breaks[2]] <- b[1]
      bij[1:breaks[2],(breaks[2]+1):breaks[3]] <- b[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- b[2]
      bij[1:breaks[4],(breaks[3]+1):breaks[4]] <- b[3]
      bij[(breaks[3]+1):breaks[4],1:breaks[4]] <- b[3]
      bij[1:breaks[5],(breaks[4]+1):breaks[5]] <- b[4]
      bij[(breaks[4]+1):breaks[5],1:breaks[5]] <- b[4]
      bij[1:breaks[6],(breaks[5]+1):breaks[6]] <- b[5]
      bij[(breaks[5]+1):breaks[6],1:breaks[6]] <- b[5]
      bij[1:breaks[7],(breaks[6]+1):breaks[7]] <- b[6]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- b[6]
    }

    if (waifw==3) {
      bij <- matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[4],1:breaks[2]] <- b[1]
      bij[1:breaks[2],1:breaks[4]] <- b[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- b[2]
      bij[(breaks[2]+1):breaks[4],(breaks[3]+1):breaks[4]] <- b[3]
      bij[(breaks[3]+1):breaks[4],(breaks[2]+1):breaks[4]] <- b[3]
      bij[1:breaks[5],(breaks[4]+1):breaks[5]] <- b[4]
      bij[(breaks[4]+1):breaks[5],1:breaks[5]] <- b[4]
      bij[1:breaks[6],(breaks[5]+1):breaks[6]] <- b[5]
      bij[(breaks[5]+1):breaks[6],1:breaks[6]] <- b[5]
      bij[1:breaks[7],(breaks[6]+1):breaks[7]] <- b[6]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- b[6]
    }

    if (waifw==4) {
      bij <- matrix(NA,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[7]] <- b[1]
      bij[(breaks[2]+1):breaks[3],1:breaks[7]] <- b[2]
      bij[(breaks[3]+1):breaks[4],1:breaks[7]] <- b[3]
      bij[(breaks[4]+1):breaks[5],1:breaks[7]] <- b[4]
      bij[(breaks[5]+1):breaks[6],1:breaks[7]] <- b[5]
      bij[(breaks[6]+1):breaks[7],1:breaks[7]] <- b[6]
    }

    if (waifw==5) {
      bij <- matrix(b[6],ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- b[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- b[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- b[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- b[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- b[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- b[5]
    }

    if (waifw==6) {
      bij <- matrix(0,ncol=Lmax,nrow=Lmax)
      bij[1:breaks[2],1:breaks[2]] <- b[1]
      bij[(breaks[2]+1):breaks[3],(breaks[2]+1):breaks[3]] <- b[2]
      bij[(breaks[3]+1):breaks[4],(breaks[3]+1):breaks[4]] <- b[3]
      bij[(breaks[4]+1):breaks[5],(breaks[4]+1):breaks[5]] <- b[4]
      bij[(breaks[5]+1):breaks[6],(breaks[5]+1):breaks[6]] <- b[5]
      bij[(breaks[6]+1):breaks[7],(breaks[6]+1):breaks[7]] <- b[6]
    }

    foiiprev <- base::rep(pwcrate(y.var=y,x.var=a,n.var=n,breaks=breaks)$ratevec, base::c(breaks[2], base::diff(breaks[-1])))
    tol <- 1
    it <- 0

    while ((tol>1e-15)&(it<5000)) {
      foii <- (N/L)*D*bij%*%(base::as.matrix(foiiprev/(foiiprev+muy))*matrix(base::c(1-exp(-(1-A)*(foiiprev[1]+muy[1])),base::exp(-(1-A)*(foiiprev[1]+muy[1])-c(0,base::cumsum(foiiprev[-1]+muy[-1])[1:(Lmax-2)]))-base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(base::cumsum(foiiprev[-1]+muy[-1])))),ncol=1))
      foii <- base::apply(base::cbind(0,foii),1,max)
      foii <- base::apply(base::cbind(1,foii),1,min)
      tol <- base::sum((foii-foiiprev)^2)
      it <- it+1
      foiiprev <- foii
    }

    if (plots=="TRUE") {
      par(mfrow=base::c(1,1))
      par(mar=base::c(5,4,4,4)+0.3)
      plot(base::c(A,1:base::max(base::floor(a))),1-base::exp(base::c(0,-(1-A)*foii[1],-(1-A)*foii[1]-base::cumsum(foii[-1])[1:(base::max(base::floor(a))-1)])),type="l",xlab="age",ylab="prevalence",ylim=base::c(0,1),xlim=base::c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),1-base::exp(-(1-A)*foii[1]-base::cumsum(foii[-1])[base::max(base::floor(a)):(Lmax-2)]),lty=2,lwd=2)
      htab <- base::table(base::floor(a),y)
      points(base::c(A,base::sort(base::unique(base::floor(a)))[-1]),htab[,2]/(htab[,1]+htab[,2]),cex=0.02*(htab[,1]+htab[,2]),lwd=1.1)
      par(new=TRUE)
      plot(base::c(A,1:base::max(base::floor(a))),foii[1:(base::max(base::floor(a))+1)],type="l",axes=FALSE,bty="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),foii[(base::max(floor(a))+2):Lmax],lty=2,lwd=2)
      axis(4,at=base::pretty(base::range(foii)))
      #contour(c(0:(Lmax-1)),c(0:(Lmax-1)),bij,xlab="age of susceptible",ylab="age of infectious")
    }

    prev <- base::rep(NA, base::length(a))
    ll <- rep(NA, base::length(a))

    for (i in 1:base::length(a)) {
      prev[i] <- (1-base::exp(base::c(-(a[i]-A)*foii[1],-(1-A)*foii[1]-base::cumsum(base::c(0,foii[-1]))-(foii[-1])[base::floor(a[i])]*(a[i]-base::floor(a[i])))))[base::floor(a[i])+1]
      ll[i] <- y[i]*base::log(prev[i]+1e-8)+(n[i]-y[i])*base::log(1-prev[i]+1e-8)
    }

    R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
    Mij <- base::diag(base::c(My))
    R0vec <- base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values

    return(list(ll=-2*base::sum(ll),eivalues=R0vec,prev=prev,bij=bij))
  }

  waifwproc.fitter <- function(b) { return(waifwproc(a=a,y=y,n=n,b=b,waifw=waifw.choice,muy=muy,breaks=breaks,Lmax=Lmax,N=N,D=D)$ll) }
  maxna <- function(x) { return(base::max(x,na.rm=TRUE)) }

  if (is.null(startpar)) {
    if (waifw.choice==1) { startpar<-base::apply(base::cbind(waifw6.fit$w1$beta,base::rep(0,6)),1,maxna) }
    if (waifw.choice==2) { startpar<-base::apply(base::cbind(waifw6.fit$w2$beta,base::rep(0,6)),1,maxna) }
    if (waifw.choice==3) { startpar<-base::apply(base::cbind(waifw6.fit$w3$beta,base::rep(0,6)),1,maxna) }
    if (waifw.choice==4) { startpar<-base::apply(base::cbind(waifw6.fit$w4$beta,base::rep(0,6)),1,maxna) }
    if (waifw.choice==5) { startpar<-base::apply(base::cbind(waifw6.fit$w5$beta,base::rep(0,6)),1,maxna) }
    if (waifw.choice==6) { startpar<-base::apply(base::cbind(waifw6.fit$w6$beta,base::rep(0,6)),1,maxna) }
  }

  waifw.result <- stats::nlm(waifwproc.fitter, base::sqrt(startpar), print.level = 0)
  result.global <- waifwproc(a=a,y=y,b=waifw.result$estimate,waifw=waifw.choice,muy=muy,breaks=breaks,Lmax=Lmax,N=N,D=D)

  return(base::list(deviance=result.global$ll,aic=result.global$ll+6*2,bic=result.global$ll+6*base::log(base::length(y)),bij=result.global$bij,R0=base::max(as.double(result.global$eivalues)),N=N,D=D,L=L))
}

contact.fitter <- function(a,y,rij,muy,N,D,Lmax,A,plots="TRUE",startpar) {
  L <- Lmax*base::mean(base::exp(-base::cumsum(muy)))
  My <- base::exp(-base::cumsum(muy))

  qproc <- function(a,y,qpar,rij,Lmax,N,D,plots="TRUE") {
    if (Lmax>100) { return("Please specify Lmax<100") }

    bij <- 365*qpar*(rij)[1:Lmax,1:Lmax]
    foiiprev <- base::rep(0.01,Lmax)
    muy <- muy[1:Lmax]
    tol <- 1
    it <- 0

    while ((tol>1e-10)&(it<2000)) {
      foii <- (N/L)*D*bij%*%(base::as.matrix(foiiprev/(foiiprev+muy))*base::matrix(base::c(1-base::exp(-(1-A)*(foiiprev[1]+muy[1])),base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(0,base::cumsum(foiiprev[-1]+muy[-1])[1:(Lmax-2)]))-base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(base::cumsum(foiiprev[-1]+muy[-1])))),ncol=1))
      foii <- base::apply(base::cbind(0,foii),1,max)
      foii <- base::apply(cbind(1,foii),1,min)
      tol <- base::sum((foii-foiiprev)^2)
      it <- it+1
      foiiprev <- foii
    }

    if (plots=="TRUE") {
      par(mfrow=c(1,1))
      par(mar=c(5,4,4,4)+0.3)
      plot(c(A,1:base::max(base::floor(a))),1-base::exp(base::c(0,-(1-A)*foii[1],-(1-A)*foii[1]-base::cumsum(foii[-1])[1:(base::max(base::floor(a))-1)])),type="l",xlab="age",ylab="prevalence",ylim=c(0,1),xlim=c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),1-base::exp(-(1-A)*foii[1]-base::cumsum(foii[-1])[base::max(base::floor(a)):(Lmax-2)]),lty=2,lwd=2)
      htab <- table(base::floor(a),y)
      points(base::c(A,base::sort(base::unique(base::floor(a)))[-1]),htab[,2]/(htab[,1]+htab[,2]),cex=0.02*(htab[,1]+htab[,2]),lwd=1.1)
      par(new=TRUE)
      plot(base::c(A,1:base::max(base::floor(a))),foii[1:(base::max(base::floor(a))+1)],type="l",axes=FALSE,bty="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),foii[(base::max(base::floor(a))+2):Lmax],lty=2,lwd=2)
      axis(4,at=base::pretty(base::range(foii)))
    }

    prev <- base::rep(NA,base::length(a))
    ll <- base::rep(NA,base::length(a))

    for (i in 1:length(a)) {
      prev[i] <- (1-base::exp(base::c(-(a[i]-A)*foii[1],-(1-A)*foii[1]-base::cumsum(base::c(0,foii[-1]))-(foii[-1])[base::floor(a[i])]*(a[i]-base::floor(a[i])))))[base::floor(a[i])+1]
      ll[i] <- y[i]*base::log(prev[i]+1e-8)+(1-y[i])*base::log(1-prev[i]+1e-8)
    }

    R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
    Mij <- base::diag(base::c(My[1:Lmax]))
    R0vec <- base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values

    return(list(ll=-2*sum(ll),eivalues=R0vec,prev=prev,bij=bij))
  }

  qproc.fitter <- function(qpar) { return(qproc(a,y,qpar,rij,Lmax,N,D,plots="TRUE")$ll) }

  q.result <- stats::nlm(qproc.fitter,startpar)
  result.global <- qproc(a=a,y=y,q=q.result$estimate,rij=rij,Lmax=Lmax,N=N,D=D)

  return(base::list(qhat=q.result$estimate,deviance=q.result$minimum,aic=q.result$minimum+2,bic=q.result$minimum+base::log(base::length(y)),bij=result.global$bij,R0=max(as.double(result.global$eivalues))))
}

contact.fitter.location <- function(a, y, rij, rij1, rij2, rij3, rij4, rij5, rij6, muy, N, D, Lmax, A, plots="TRUE", startpar) {
  no.rij <- base::max((rij1!=0))+base::max((rij2!=0))+base::max((rij3!=0))+base::max((rij4!=0))+base::max((rij5!=0))+base::max((rij6!=0))
  L <- Lmax*base::mean(base::exp(-base::cumsum(muy)))
  My <- base::exp(-base::cumsum(muy))

  qproc <- function(a,y,qpar,rij,Lmax,N,D,plots="TRUE") {
    qpar<-qpar^2

    if (Lmax>100) { return("Please specify Lmax<100") }

    bij <- 365*(qpar[1]*(rij1)[1:Lmax,1:Lmax]+qpar[2]*(rij2)[1:Lmax,1:Lmax]+qpar[3]*(rij3)[1:Lmax,1:Lmax]+qpar[4]*(rij4)[1:Lmax,1:Lmax]+qpar[5]*(rij5)[1:Lmax,1:Lmax]+qpar[6]*(rij6)[1:Lmax,1:Lmax])
    foiiprev <- base::rep(0.01,Lmax)
    muy <- muy[1:Lmax]
    tol <- 1
    it <- 0

    while ((tol>1e-10)&(it<2000)) {
      foii <- (N/L)*D*bij%*%(base::as.matrix(foiiprev/(foiiprev+muy))*base::matrix(base::c(1-base::exp(-(1-A)*(foiiprev[1]+muy[1])),base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(0,base::cumsum(foiiprev[-1]+muy[-1])[1:(Lmax-2)]))-base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(base::cumsum(foiiprev[-1]+muy[-1])))),ncol=1))
      foii <- base::apply(base::cbind(0,foii),1,max)
      foii <- base::apply(base::cbind(1,foii),1,min)
      tol <- base::sum((foii-foiiprev)^2)
      it <- it+1
      foiiprev <- foii
    }

    if (plots=="TRUE") {
      par(mfrow=base::c(1,1))
      par(mar=base::c(5,4,4,4)+0.3)
      plot(base::c(A,1:base::max(base::floor(a))),1-base::exp(base::c(0,-(1-A)*foii[1],-(1-A)*foii[1]-base::cumsum(foii[-1])[1:(base::max(base::floor(a))-1)])),type="l",xlab="age",ylab="prevalence",ylim=base::c(0,1),xlim=base::c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),1-base::exp(-(1-A)*foii[1]-base::cumsum(foii[-1])[base::max(base::floor(a)):(Lmax-2)]),lty=2,lwd=2)
      htab <- base::table(base::floor(a),y)
      points(base::c(A,base::sort(base::unique(base::floor(a)))[-1]),htab[,2]/(htab[,1]+htab[,2]),cex=0.02*(htab[,1]+htab[,2]),lwd=1.1)
      par(new=TRUE)
      plot(base::c(A,1:base::max(base::floor(a))),foii[1:(base::max(base::floor(a))+1)],type="l",axes=FALSE,bty="n",xlab="",ylab="",ylim=base::c(0,1),xlim=base::c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),foii[(base::max(base::floor(a))+2):Lmax],lty=2,lwd=2)
      axis(4,at=base::pretty(base::range(foii)))
    }

    prev <- base::rep(NA,base::length(a))
    ll <- base::rep(NA,base::length(a))

    for (i in 1:base::length(a)) {
      prev[i] <- (1-base::exp(base::c(-(a[i]-A)*foii[1],-(1-A)*foii[1]-base::cumsum(base::c(0,foii[-1]))-(foii[-1])[base::floor(a[i])]*(a[i]-base::floor(a[i])))))[base::floor(a[i])+1]
      ll[i] <- y[i]*base::log(prev[i]+1e-8)+(1-y[i])*base::log(1-prev[i]+1e-8)
    }

    R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
    Mij <- base::diag(base::c(My[1:Lmax]))
    R0vec <- base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values

    return(list(ll=-2*base::sum(ll),eivalues=R0vec,prev=prev,bij=bij))
  }

  qproc.fitter <- function(qpar) { return(qproc(a,y,qpar,rij,Lmax,N,D,plots="TRUE")$ll) }

  q.result <- stats::nlm(qproc.fitter,base::sqrt(startpar),hessian=TRUE)
  result.global <- qproc(a=a,y=y,q=q.result$estimate,rij=rij,Lmax=Lmax,N=N,D=D)

  return(base::list(qhat=q.result$estimate^2,qhess=q.result$hessian,deviance=q.result$minimum,aic=q.result$minimum+no.rij*2,bic=q.result$minimum+no.rij*base::log(length(y)),bij=result.global$bij,R0=base::max(base::as.double(result.global$eivalues))))
}

contact.fitter.loglinear <- function(a, y, rij, int=FALSE, muy, N, D, Lmax, A, plots="TRUE", startpar) {
  L <- Lmax*base::mean(base::exp(-base::cumsum(muy)))
  My <- base::exp(-base::cumsum(muy))

  qproc <- function(a,y,qpar,rij,Lmax,N,D,plots="TRUE") {
    if (int==FALSE) { qpar[3]=0 }
    if (Lmax>100) { return("Please specify Lmax<100") }

    q.f <- function(x,y) { base::exp(qpar[1]+qpar[2]*x+qpar[2]*y+qpar[3]*y*x) }

    qij <- base::outer(base::c(1:Lmax),c(1:Lmax),q.f)
    bij <- 365*qij*(rij)[1:Lmax,1:Lmax]
    foiiprev <- base::rep(0.01,Lmax)
    muy <- muy[1:Lmax]
    tol <-1
    it <-0

    while ((tol>1e-10)&(it<2000)) {
      foii <- (N/L)*D*bij%*%(base::as.matrix(foiiprev/(foiiprev+muy))*base::matrix(base::c(1-base::exp(-(1-A)*(foiiprev[1]+muy[1])),base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(0,base::cumsum(foiiprev[-1]+muy[-1])[1:(Lmax-2)]))-base::exp(-(1-A)*(foiiprev[1]+muy[1])-base::c(base::cumsum(foiiprev[-1]+muy[-1])))),ncol=1))
      foii <- base::apply(cbind(0,foii),1,max)
      foii <- base::apply(cbind(1,foii),1,min)
      tol <- base::sum((foii-foiiprev)^2)
      it <- it+1
      foiiprev <- foii
    }

    if (plots=="TRUE") {
      par(mfrow=c(1,1))
      par(mar=base::c(5,4,4,4)+0.3)
      plot(base::c(A,1:base::max(base::floor(a))),1-base::exp(base::c(0,-(1-A)*foii[1],-(1-A)*foii[1]-base::cumsum(foii[-1])[1:(base::max(base::floor(a))-1)])),type="l",xlab="age",ylab="prevalence",ylim=base::c(0,1),xlim=base::c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),1-base::exp(-(1-A)*foii[1]-base::cumsum(foii[-1])[base::max(base::floor(a)):(Lmax-2)]),lty=2,lwd=2)
      htab <- base::table(base::floor(a),y)
      points(base::c(A,base::sort(base::unique(base::floor(a)))[-1]),htab[,2]/(htab[,1]+htab[,2]),cex=0.02*(htab[,1]+htab[,2]),lwd=1.1)
      par(new=TRUE)
      plot(c(A,1:base::max(base::floor(a))),foii[1:(base::max(base::floor(a))+1)],type="l",axes=FALSE,bty="n",xlab="",ylab="",ylim=base::c(0,1),xlim=base::c(0,80),lwd=2)
      lines((base::max(base::floor(a))+1):(Lmax-1),foii[(base::max(base::floor(a))+2):Lmax],lty=2,lwd=2)
      axis(4,at=base::pretty(base::range(foii)))
    }

    prev <- base::rep(NA,base::length(a))
    ll <- base::rep(NA,base::length(a))

    for (i in 1:base::length(a)) {
      prev[i]<-(1-base::exp(base::c(-(a[i]-A)*foii[1],-(1-A)*foii[1]-base::cumsum(base::c(0,foii[-1]))-(foii[-1])[base::floor(a[i])]*(a[i]-base::floor(a[i])))))[base::floor(a[i])+1]
      ll[i]<-y[i]*base::log(prev[i]+1e-8)+(1-y[i])*base::log(1-prev[i]+1e-8)
    }

    R0ij <- (N/L)*D*bij[1:Lmax,1:Lmax]
    Mij <- base::diag(base::c(My[1:Lmax]))

    if (base::sum(base::is.na(Mij%*%R0ij))==0) { R0vec<- base::eigen(Mij%*%R0ij,symmetric=FALSE,only.values=TRUE,EISPACK=FALSE)$values }

    return(base::list(ll=-2*base::sum(ll),eivalues=R0vec,prev=prev,bij=bij))
  }

  qproc.fitter <- function(qpar){ return(qproc(a=a,y=y,qpar=qpar,rij=rij,Lmax=Lmax,N=N,D=D,plots="TRUE")$ll) }

  q.result <- stats::nlm(qproc.fitter,startpar,hessian=TRUE)
  result.global <- qproc(a=a,y=y,q=q.result$estimate,rij=rij,Lmax=Lmax,N=N,D=D)

  return(base::list(qhat=q.result$estimate[q.result$estimate!=0],qhess=q.result$hessian[q.result$estimate!=0,q.result$estimate!=0],deviance=q.result$minimum,aic=q.result$minimum+base::sum(q.result$estimate!=0)*2,bic=q.result$minimum+base::sum(q.result$estimate!=0)*base::log(base::length(y)),bij=result.global$bij,L=L,D=D,N=N,R0=base::max(base::as.double(result.global$eivalues))))
}
