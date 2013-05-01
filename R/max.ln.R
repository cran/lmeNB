max.ln <-
function(x=1, y1=0, maxY2=3, u1=1.5, u2=c(1.5,1.5,1.5), a_inv=0.5, 
       #see max.fun  
       sh=0.5, sc=2 #sh=mean, sc=sd of the lognormal dist'n
    )
{   
    P=x/(x+a_inv) 
    #pr(Y1+= y1+ |g) 
    tem=P^y1*(1-P)^u1*dlnorm(x, meanlog=sh, sdlog=sc)
    for (i in 1:length(u2))
    { tem1=pnbinom(maxY2-1, prob=1-P, size=u2[i])
      tem=tem*tem1
    }
    return(tem)
}
