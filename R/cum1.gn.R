cum1.gn <-
function(x=2, a_inv=0.5, y1=2, y2=2, u1=3, u2=3, 
                                             #see cum.fun
                 sh=0.45, sc=2,             #gamma parameters 
                 u.n=3, s.n=0.5, p.mx=0.05  # normal parameters and proportion, see getSH.gn
                )
{   p=x/(x+a_inv)  
    dg=p.mx*dnorm(x, mean=u.n, sd=s.n)+(1-p.mx)*dgamma(x, shape=sh, scale=sc)
    tem=p^y1*(1-p)^u1*dg
    tem=tem*pnbinom(y2-1, prob=1-p, size=u2)
    return(tem)
}




cum1.uf <-
function(x=0, a_inv=0.5, y1=2, y2=2, u1=3, u2=3) #see cum.fun
{   x1=exp(x)
    p=x1/(x1+a_inv)  
    tem=p^y1*(1-p)^u1
    tem=tem*pnbinom(y2-1, prob=1-p, size=u2)
    return(tem)
}
