dmarg.fun <-
function(y=1, u=1.5, a=0.5, th=3, dist="G")
{ uVa=u/a
  # "G": G ~ gamma
  if (dist=="G") tem1=integrate(dint.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/th, sc=th, ysum=y, usum=uVa)
  # "N": G ~ Log normal
  if (dist=="N")  {
   tem=log(th+1)
   u.ln=-tem/2
   s.ln=sqrt(tem)
   tem1=integrate(dint.ln.fun, lower=0, upper=Inf, a_inv=1/a, sh=u.ln, sc=s.ln, ysum=y, usum=uVa)}

  val=tem1$v
  return(val)
}

dmarg.gauss.fun <-
function(y=1, u=1.5, a=0.5, th=3, dist="G", gau.n=32)
{   ainv=1/a 
    uVa=u/a

    if (dist=="G")  
    { sc=th 
      sh=1/th
      GAU <- gauss.quad.prob(n=gau.n, dist="gamma", alpha=sh, beta=sc)
    }
    if (dist=="N")  
    { tem=log(th+1)
      u.ln=-tem/2
      s.ln=sqrt(tem)
      GAU <- gauss.quad.prob(n=gau.n,dist="normal", mu=u.ln, sigma=s.ln)
      GAU$nodes=exp(GAU$nodes)   
    }

    p=GAU$nodes/(GAU$nodes+ainv)
   
    val=sum(dnbinom(y, size=uVa, prob=1-p)*GAU$weights)
    return(val)
}

