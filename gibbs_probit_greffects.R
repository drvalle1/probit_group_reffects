rm(list=ls())
library('mvtnorm')
set.seed(35)

setwd('U:\\GIT_models\\probit_group_reffects')
source('aux probit greffects.R')

setwd('U:\\GIT_models\\probit_group_reffects\\fake data')
y=data.matrix(read.csv('fake data y.csv',as.is=T))
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
xtx=t(xmat)%*%xmat

#basic settings
nobs=nrow(y)
nspp=ncol(y)
ngr=10
nparam=ncol(xmat)

#priors
prior.a=10; prior.b=1 #mean of precision 10, means that tau2=0.1
gamma1=0.1

#initial values
mu=matrix(0,nparam,ngr)
tau2=matrix(1,nparam,ngr)
w=sample(1:ngr,size=nspp,replace=T)
betas=matrix(0,nparam,nspp)
theta=rep(1/ngr,ngr)
ystar=matrix(ifelse(y==1,1,-1),nobs,nspp)

#MCMC settings
ngibbs=1000; nburn=ngibbs/2
store.tau2=store.mu=matrix(NA,ngibbs,nparam*ngr)
store.w=matrix(NA,ngibbs,nspp)
store.betas=matrix(NA,ngibbs,nparam*nspp)
store.theta=matrix(NA,ngibbs,ngr)
options(warn=2)
for (i in 1:ngibbs){
  print(i)
  print(max(w))
  
  mu=sample.mu(w=w,betas=betas,tau2=tau2,nparam=nparam,ngr=ngr)
  tau2=sample.tau2(betas=betas,mu=mu,prior.a=prior.a,prior.b=prior.b,w=w)
  w=sample.w(tau2=tau2,betas=betas,ltheta=log(theta),w=w,ngr=ngr,mu=mu)
  betas=sample.betas(mu=mu,tau2=tau2,xmat=xmat,
                     ystar=ystar,nspp=nspp,nparam=nparam,
                     w=w,xtx=xtx)
  theta=sample.theta(gamma1=gamma1,w=w,ngr=ngr)
  # theta=rep(1/ngr,ngr)
  ystar=sample.ystar(y=y,xmat=xmat,betas=betas,nobs=nobs,nspp=nspp)
  
  #re-order w from time to time
  if (i<nburn & i%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    mu=mu[,ind]
    tau2=tau2[,ind]
    wnew=w
    for (j in 1:ngr){
      cond=w==ind[j]
      wnew[cond]=j
    }
    w=wnew
  }

  #store results
  store.tau2[i,]=tau2
  store.mu[i,]=mu
  store.w[i,]=w
  store.betas[i,]=betas
  store.theta[i,]=theta
}
