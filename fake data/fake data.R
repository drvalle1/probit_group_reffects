rm(list=ls())
set.seed(32)

#basic settings
nobs=1000
nspp=150
ngr=5
nparam=4

#parameters
seq1=seq(from=-2,to=2,by=1)
mu.true=mu=matrix(sample(seq1,size=nparam*ngr,replace=T),nparam,ngr)
tau.true=tau=matrix(runif(nparam*ngr,min=0,max=0.2),nparam,ngr)
w.true=w=sample(1:ngr,size=nspp,replace=T)

#get betas
betas=matrix(NA,nparam,nspp)
for (i in 1:nspp){
  mu1=mu[,w[i]]
  tau1=tau[,w[i]]
  betas[,i]=rnorm(nparam,mean=mu1,sd=tau1)
}
betas.true=betas

mu.estim=numeric()
for (i in 1:ngr){
  tmp=rowMeans(betas[,w==i])  
  mu.estim=cbind(mu.estim,tmp)
}
plot(mu.estim,mu.true)

#visualize these results
rango=range(betas)
par(mfrow=c(3,1),mar=rep(1,4))
for (i in 1:ngr) image(betas[,w==i],zlim=rango)

#get covariates
xmat=matrix(rnorm(nobs*nparam),nobs,nparam)
xmat[,1]=1

#generate observations
media=xmat%*%betas
ystar=matrix(rnorm(nobs*nspp,mean=media,sd=1),nobs,nspp)
y=matrix(ifelse(ystar>0,1,0),nobs,nspp)

#export data
setwd('U:\\GIT_models\\probit_group_reffects\\fake data')
write.csv(y,'fake data y.csv',row.names=F)
write.csv(xmat,'fake data xmat.csv',row.names=F)
