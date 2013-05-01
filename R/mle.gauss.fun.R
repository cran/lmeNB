mle.gauss.fun <-
function(dat, p.ini, IPRT=FALSE, model="G", #see mle.fun
                       gau.n=64  # number of nodes for gauss quad
                      ) 
{ 
  DT=getDT(dat)
  if (model=="G") mod="gamma"   
  if (model=="N") mod="normal"

  tt=optim(p.ini, lk.gauss.fun, hessian=T, dat=DT, Iprt=IPRT, mod=mod, gau.n=gau.n)
    
  nlk=tt$value+sum(lgamma(DT$y+1))
  if (IPRT) print(tt$hessian)
  #vcm=vcm.fun(he=tt$hessian, corr=F)
  vcm=solve(tt$hessian)
  p.est=cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est)=c("log_a", "log_th", "log_u0", DT$xnames)
  return(list(opt=tt, nlk=nlk, V=vcm, est=p.est, mod=model))
}


mle.pois.fun <-
function(dat, p.ini, IPRT=FALSE, model="G") #G = gamma 
{ 
  DT=getDT(dat)
  if (model=="G")
  { tt=optim(p.ini, lk.pois.fun, hessian=T, dat=DT, Iprt=IPRT)
  }
#  if (model=="N")
#  {   
#      tt=optim(p.ini, lk.pois.ln.fun, hessian=T, dat=DT, Iprt=IPRT, g.nodes=64)
#  }

  nlk=tt$value+sum(lgamma(DT$y+1))
  if (IPRT) print(tt$hessian)
  vcm=solve(tt$hessian)
  p.est=cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est)=c("log_th", "log_u0", DT$xnames)

  return(list(opt=tt, nlk=nlk, V=vcm, est=p.est, mod=model))
}
