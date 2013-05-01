lk.gauss.fun <-
function(para, dat, Iprt=F, gau.n=64, mod="gamma", sig=F)
{ if (Iprt) cat(para)
  a=exp(para[1])
  th1=exp(para[2])  #scale parameter for dgamma()
  ainv=1/a  
  if (mod=="gamma") shp=1/th1

  if (mod=="normal") 
  { s.ln=log(th1+1) # s^2
    u.ln=-s.ln/2    #-s^2/2 #mean for dlnorm
    s.ln=sqrt(s.ln) #sd for dlnorm()
  }

  u0=exp(para[3])
  th2=rep(u0/a, dat$totN)
  if (dat$cn>0) {
   b=para[4:(dat$cn+3)]
   tem=exp(dat$x%*%b)
   th2=tem*th2
  }

   nllk=sum( - lgamma(dat$y+th2) + lgamma(th2))
   us=tapply(th2, dat$id, sum)

   #nodes for gauss quad
   if (mod=="gamma")  
   { GAU <- gauss.quad.prob(n=gau.n,dist=mod, alpha=shp, beta=th1)}
   if (mod=="normal")  
   { GAU <- gauss.quad.prob(n=gau.n,dist=mod, mu=u.ln, sigma=s.ln)
     GAU$nodes=exp(GAU$nodes)   
   }

   p=GAU$nodes/(GAU$nodes+ainv)
   q=1-p

   lk=rep(0, dat$np)
   for (i in 1:dat$np)
   { tem = sum(GAU$weights*p^dat$ys[i]*q^us[i])

     #print(c(dat$ys[i], us[i], tem))
     nllk = nllk - log(tem)
     if (sig)  #-loglikelihood for each patient
     {  ll=dat$ind[i]:(dat$ind[i+1]-1)
        lk[i]=log(tem)+sum(lgamma(dat$y[ll]+th2[ll]) - lgamma(th2[ll])-lgamma(dat$y[ll]+1))
     }
   }
   if (Iprt) cat(" nllk=", nllk, "\n")
   if (sig) return(lk) 
   else return(nllk)
}


lk.ln.fun <-
function(para, dat, Iprt=F, tol=1e-75)
{ if (Iprt) cat(para)
  
  a=exp(para[1])
  th1=exp(para[2]) #scale
  
  s.ln=log(th1+1) #s^2
  u.ln=-s.ln/2    #-s^2/2, mean for dlnorm()
  s.ln=sqrt(s.ln) #s, sd for dlnorm()

  ainv=1/a  

  u0=exp(para[3])
  th2=rep(u0/a, dat$totN)
  if (dat$cn>0) {
   b=para[4:(dat$cn+3)]
   tem=exp(dat$x%*%b)
   th2=tem*th2
  }

   nllk=sum( - lgamma(dat$y+th2) + lgamma(th2))
   us=tapply(th2, dat$id, sum)
   
   for (i in 1:dat$np) 
   { tem=integrate(int1.ln, lower=0, upper=Inf, a_inv=ainv,
            sh=u.ln, sc=s.ln, ysum=dat$ys[i],  usum=us[i], abs.tol=tol)
     #print(c(dat$ys[i], us[i], tem$v))
     nllk = nllk - log(tem$v)
   }
   if (Iprt) cat(" nllk=", nllk, "\n")
   return(nllk)
}


lk.pois.fun <-
function(para, dat, Iprt=F, tol=1.e-75, sig=F)
{ if (Iprt) cat(para)
  th1=exp(para[1]) 
  
  u=rep(para[2], dat$totN)
  if (dat$cn>0) {
   b=para[3:(dat$cn+2)]
   u=u+dat$x%*%b
  }
 th2=exp(u)

   us=tapply(th2, dat$id, sum)
   thy= th1+dat$ys

   #log likelihood without constant
   nllk= - sum(dat$y*u) - sum(lgamma(thy) -  thy*log(th1+us)) - dat$np*(th1*para[1] - lgamma(th1))

   #sig=T -loglikelihood for each subject
   if (sig) 
   {  nllk = nllk + sum(lgamma(dat$y+1))
   }
   if (Iprt) cat(" nllk=", nllk, "\n")
  return(nllk)   
}
