compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango,col='red',lwd=2)  
}

theta.estim=store.theta[ngibbs,]
plot(theta.estim,type='h')

w.estim=store.w[ngibbs,]
aux=data.frame(w.true=w.true,w.estim=w.estim)
k=table(aux); k
seq1=c(3,1,2,5,4)
k[,seq1]

betas.estim=matrix(store.betas[ngibbs,],nparam,nspp)
compare1(estim=betas.estim,true=betas.true)

ngr=10
mu.estim=matrix(store.mu[ngibbs,],nparam,ngr)
compare1(estim=mu.estim[,seq1],true=mu.true)

tau2.estim=matrix(store.tau2[ngibbs,],nparam,ngr)
compare1(estim=tau2.estim[,seq1],true=tau.true^2)