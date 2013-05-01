pint.fun <-
function(x=2, a_inv=0.5, sh=0.5, sc=2, ysum=2, usum=3)
{   p=x/(x+a_inv) 
    tem=pnbinom(ysum, size=usum, prob=1-p)*dgamma(x, shape=sh, scale=sc)
    return(tem)
}
