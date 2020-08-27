tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
sample.mu=function(w,betas,tau2,nparam,ngr){
  mu=matrix(NA,nparam,ngr)
  for (i in 1:ngr){
    cond=w==i
    nw=sum(cond)
    invT=diag(1/tau2[,i])
    prec=(1/100)+(nw/tau2[,i])
    var1=diag(1/prec)
    if (nw==0) soma.betas=matrix(0,nparam,1)
    if (nw==1) soma.betas=betas[,cond]
    if (nw >1) soma.betas=matrix(rowSums(betas[,cond]),nparam,1)
    pmedia=invT%*%soma.betas
    mu[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
  }
  mu
}
sample.tau2=function(betas,mu,prior.a,prior.b,w){
  tau2=matrix(NA,nparam,ngr)
  for (i in 1:ngr){
    cond=w==i
    nw=sum(cond)
    a1=(nw+2*prior.a)/2

    if (nw==0) err2=0
    if (nw==1) err2=(betas[,cond]-mu[,i])^2
    if (nw >1) {
      mu1=matrix(mu[,i],nparam,nw)
      err2=rowSums((betas[,cond]-mu1)^2)
    }
    b1=(err2/2)+prior.b
    tau2[,i]=1/rgamma(nparam,a1,b1)
  }
  tau2
}
sample.w=function(tau2,betas,ltheta,w,ngr,mu){
  p1=(-1/2)*colSums(log(tau2))
  for (i in 1:nspp){
    err=(betas[,i]-mu)^2
    var.part=(-1/(2*tau2))
    p2=colSums(var.part*err)
    lprob=p1+p2+ltheta

    #do we need a new group?    
    max.w=max(w)
    if (max.w<ngr){
      lprob=lprob[1:max.w]

      #probability of new group
      ind=max.w+1
      p1=(-1/2)*sum(log(100+tau2[,ind]))
      err2=betas[,ind]^2
      var.part=-1/(2*(100+tau2[,ind]))
      p2=sum(var.part*err2)
      lprob=c(lprob,c(p1+p2+ltheta[ind]))
    }
    tmp=lprob-max(lprob)
    tmp=exp(tmp)
    prob=tmp/sum(tmp)
    ind=rmultinom(1,size=1,prob=prob)
    w[i]=which(ind==1)
  }
  w
}
sample.betas=function(mu,tau2,xmat,ystar,nspp,nparam,w,xtx){
  betas=matrix(NA,nparam,nspp)
  prec=1/tau2
  for (i in 1:nspp){
    w1=w[i]
    invT=diag(prec[,w1])
    prec1=xtx+invT
    var1=solve(prec1)
    pmedia=t(xmat)%*%ystar[,i]+invT%*%mu[,w1]
    betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
  }
  betas
}
sample.theta=function(gamma1,w,ngr){
  nk=rep(0,ngr)
  tmp=table(w)
  nk[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,ngr)
  aux=1
  for(i in 1:(ngr-1)){
    nk.maior=nk[(i+1):ngr]
    v[i]=rbeta(1,nk[i]+1,sum(nk.maior)+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  theta[ngr]=aux
  theta
}
sample.ystar=function(y,xmat,betas,nobs,nspp){
  media=xmat%*%betas
  lo.mat=matrix(ifelse(y==1,0,-1000),nobs,nspp)
  hi.mat=matrix(ifelse(y==0,0, 1000),nobs,nspp)
  ystar=tnorm(n=nobs*nspp,lo=lo.mat,hi=hi.mat,mu=media,sig=1)
  matrix(ystar,nobs,nspp)
}