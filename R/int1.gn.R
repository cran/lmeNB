int1.gn <-
function(x=2, a_inv=0.5, sh=0.45, sc=2, u.n=3, s.n=0.5, p.mx=0.05, ysum=2, usum=3)
#see cum1.gn
{   p=x/(x+a_inv)  
    dg=p.mx*dnorm(x, mean=u.n, sd=s.n)+(1-p.mx)*dgamma(x, shape=sh, scale=sc)
    tem=p^ysum*(1-p)^usum*dg
    return(tem)
}


int1.ln <-
function(x=2, 
                 a_inv=0.5, #see int.fun
                 sh=0.5, #mean for dlnorm()
                 sc=2,   #sd for dlnorm()
                 ysum=2, usum=3 #see int.fun
                )
{   p=x/(x+a_inv)  
    tem=p^ysum*(1-p)^usum*dlnorm(x, meanlog=sh, sdlog=sc)
    return(tem)
}
