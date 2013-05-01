Pmax.non <-
function(Y1=0,               #sum(ypre); ypre+ 
                  Y2=1,               #max(ynew)
                  u1=4.5,             #Ex(Ypre+)
                  u2=c(1.5,1.5, 1.5), #Ex(Ynew), vector
                  a=0.5, 
                  gi=NULL)            #a vector (todo) allow a freq table [freq, valu]
{
    pb=1/(gi*a+1)

    p=1
    if (Y2>0)
    { j2=j1=dnbinom(Y1, prob=pb, size=u1/a)
      for (i in 1:length(u2))
      { j2=pnbinom(Y2-1, prob=pb, size=u2[i]/a)*j2 } #pr(Ynew[i]<Y2)
      p=1-sum(j2)/sum(j1)
    }
    return(p)
}
