Psum.non <-
function(Y1=0, Y2=1,     # Y1=sum(y.pre), Y2=sum(y.new)
                  u1=1.5, u2=1.5, #u1 = Ex(Y1); u2= Ex(Y2)
                  a=0.5, 
                  gi=NULL         #a freqency table = [freq, g.value]
                 )
{
    pb=1/(gi[,2]*a+1)
    
    p=1
    if (Y2>0)
    { j1=dnbinom(Y1, prob=pb, size=u1/a)*gi[,1]
      j2=pnbinom(Y2-1, prob=pb, size=u2/a)*j1
      p=1-sum(j2)/sum(j1)
    }
    return(p)
}
