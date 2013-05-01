nr.fun <-
function(th)
{ a=log(2)+log(th+1) 
  b=exp(a)
 #need solve f1=0 and f2=0 (see nr.pdf)
 f1=exp(a)-exp(-b)-(a+b)
 f2=exp(2*a)-exp(-2*b)-2*(a+b)*(th+1)
 
 ch=max(abs(f1), abs(f2))
 while (ch>0.0001)
{
 #print(c(a,b, f1,f2) )
 df11=exp(a)-1
 df12=exp(-b)-1

 df21=2*exp(2*a)-2*th-2
 df22=2*exp(-2*b)-2*th-2

 m=matrix(c(df11,df21, df12, df22),2,2)
 v=-c(f1,f2)

 tem=solve(m,v)

 a=a+tem[1]
 b=b+tem[2]
 f1=exp(a)-exp(-b)-(a+b)
 f2=exp(2*a)-exp(-2*b)-2*(a+b)*(th+1)
 
 ch=max(abs(f1), abs(f2))
}

 #print(c(f1,f2))
 #print(c(a,b))
 return(c(a,b))
}
