                                        #required library(s): numDeriv
                                        #other files required: mle-nb-mix.R

                                        # Functions to fit the Negative binomial mixture model with an AR(1) dependence structure
                                        # current options for distribution of random effects (RE): gamma, lognormal, non-parametric


##===========MLE=====================
mle.ar1.fun <- function(
                        formula,
                        ## an object of class "formula"
                        ## (or one that can be coerced to that class):
                        ## a symbolic description of the model to be fitted.
                        data,
                        ## a data frame, list or environment (or object coercible
                        ## by as.data.frame to a data frame)
                        ## containing the variables in the model.
                        ID,
                        ## a vector of length n*ni containing patient IDs of observations in data
                        Vcode,
                        ## scan number need to be integers with increment of one, i.e., -1, 0, 1, 2,..
                        p.ini=NULL,     ## initial values for the parameters
                        ## c(log(a), log(th), lgt(dt), b0, b1, ...)
                        IPRT=FALSE,    ## FALSE control: T = print iterations
                        model="G",    # dist'n for RE G= gamma, N = lognormal
                        i.tol=1.e-75, # tolerance; for integration (ar1.lk)
                        o.tol=1.e-3  # tolerance: for optim 
                        ## enter the number of scans when there is no missing data 
                        ) 
{
  dat <- formulaToDat(formula=formula,data=data,ID=ID)
  ## dat = (ID, Y, x1, x2, ...) numeric matrix
  DT=getDT(dat) ## list

  DT$dif=c(0,diff(Vcode))   
  DT$dif[DT$ind]=0
  ## DT$ind contains the locations of the 1st repeated measures for each patient
  ## DT$diff  = 0 if its the first repeated measure and = 1 ifelse
  DT$dif=DT$dif[1:DT$totN]    # scan lag
  
  ##FixN=max(DT$ni)

  if (is.null(p.ini))
    {
      p.ini=rep(0, 4+DT$cn)
      p.ini[4]=mean(DT$y)   
    }
  
  if (IPRT)
    cat("\n\n estimates: log(a),log(theta),lgt(d),b0,b1,... and the negative of the log-likelihood")
  
  tt <- optim(p.ini,  ## c(log(a), log(th), lgt(dt), b0, b1, ...)
              ar1.lk, ##ar1.lk likelihood function for gamma/log-normal RE model 
              hessian=TRUE,  
              control=list(reltol=o.tol),
              dat=DT,#input data (output from getDT; lag)
              Iprt=IPRT,
              tol=i.tol,
              ##FixN=FixN, # FixN: if uij = uj
              dist=model ## Notice that the model option is passed to ar1.lk
    )
  
  nlk <- tt$value #neg likelihood

  vcm <- solve(tt$hessian)
  if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c("log_a", "log_th", "logit_dt","(Intercept)", DT$xnames)
  p.est <- cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est) <- c("log_a", "log_th", "logit_dt", "(Intercept)", DT$xnames)
  re <- list(opt=tt, nlk=nlk, V=vcm, est=p.est, mod=model,##idat=data.frame(dat),
              Vcode=Vcode,cor="ar1",formula=formula)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}



ar1.lk <- function(para,
                   ## c(log(a), log(th), lgt(dt), b0, b1, ...)
                   dat, #input data (output from getDT; lag)
                   Iprt=TRUE,      #print control
                   tol=1.e-75,  # tolerance for integration 
                   sig=FALSE,       # if T, the compute full likelihood
                   ##FixN,     # FixN: if uij = uj
                   dist = "G"   #dist'n of the RE; G= gamma, N = lognormal
                   ) 
{
  if(Iprt) cat("\n",para," ")
  ainv=exp(-para[1])    ## ainv=1/a  
  th1=exp(para[2])      ## scalar of gamma /var(G_i) of lognormal
  if (dist=="G") shp=1/th1           ## gamma-shape

  ## When G_i ~ LN, var(G_i)=theta
  if (dist=="N") 
    { s.ln=log(th1+1) ## sigma^2 = log(th1+1) of log-normal
      u.ln=-s.ln/2    ## mu = -log(th1+1)/2 of log-normal
      s.ln=sqrt(s.ln) # sigma of log-normal
    }

  dt=ilgt(para[3])    #inverse logit = delta
  
  tem=rep(0, dat$totN) #tem = zero vector of length (s*sn) if there is no covariate 
  ## If the number of covariates is greater than 1
  if (dat$cn>0) {
    b=para[5:(dat$cn+4)] ## = [b1 b2 ...]
    tem=dat$x%*%b ## tem = b1*x1 + b2*x2 + ... (dat$x is (n*sn) by # covariates)
  }
                                        #r[i,j] = u[i,j]/alpha =exp(- log(alpha)+b0 + b1*x1 + ... )
  th2=exp(tem+ ## b1*x1+b2*x2+...
    para[4]- ## b0
    para[1] ## log a
    ) ## th2 is a vector of length (n*sn) 
  
  Dl=dt^dat$dif ## Dl = 1 for the first repeated measure of each patient and  Dl = d ifelse

  szm=c(0, th2[-dat$totN]) #r[i, j-1]
                                        #cat("Dl=", length(Dl), "szm=", length(szm), "th2=", length(th2),"\n")
  
  U = Dl*szm  #d*r[i,j-1] 
  V = szm-U   #(1-d)*r[i, j-1]
  sz2 = th2-U #r[i, j] - d*r[i, j-1]
  tem = dat$ind[1:dat$np]
  sz2[tem] = th2[tem] #set sz2 = r[i, j] for the first scan
                                      
  if (any(sz2<=0)) return(1.e15)

  nllk=0               #total likelihood
  lki=rep(0, dat$np)   #likelihood for each individual (when sig=T)

  llk0=0               #likelihood for an observation with all 0 counts
  for (i in 1:dat$np) ## i in 1, ..., # patients
    { ll=dat$ind[i]:(dat$ind[i+1]-1)
      ## ll is a vector, containing the locations of the repeated measures
      ## corresponding to i^th obs.

      ##observations with all 0 counts
      ##if (dat$ni[i]==FixN & dat$ys[i]==0 & llk0<0) 
      ##{ lki[i]=llk0      
      ##}
      ##else
      ##{
      if (dist=="G") tem <- integrate(ar1.intg, lower=0, upper=Inf, abs.tol=tol,
            a_inv=ainv, sh=shp, sc=th1, 
            y=dat$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
      
      if (dist=="N") tem <- integrate(ar1.ln.intg, lower=0, upper=Inf, abs.tol=tol, 
            a_inv=ainv, mu=u.ln, sig=s.ln, #check
            y=dat$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
                                        #print(c(dat$ys[i], us[i], tem$v))
      lki[i]=log(tem$value)
      ##if (dat$ni[i]==FixN&dat$ys[i]==0) llk0=lki[i]
      ##}
      nllk=nllk-lki[i]
    }
  if (Iprt) cat(" nllk=", nllk, "\n")

  if (sig) return(lki) 
  else return(nllk)
}

##integradient oc to prob(Y|g)*f(g); G~gamma
ar1.intg <- function(x=2,       #G=x
                     a_inv=0.5, #1/alpha 
                     sh=0.5, sc=2, # gamma parameters
                     y=1:3,        #count 
                     u=c(0,1,1),   # d*r[i,j-1] 
                     v=c(0,2,2),   # (1-d)*r[i,j-1]
                     s2=c(3,2,2)   # r[i,j]- d*r[i, j-1]
                     )
  ##y,u,v, s2 need to be of the same length
  
{ 
  Pr=a_inv/(x+a_inv)
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem=c(tem,ar1.fun(y, u, v, s2, pr))
    }
  res=tem*dgamma(x, shape=sh, scale=sc)
  return(res)
}

                                        #example: integrate with gamma density
                                        #integrate(ar1.intg, lower=0, upper=Inf, a_inv=0.5, sh=1/5, sc=5, y=1:3)

##integradient oc to prob(Y|g)*f(g); G~lognormal
ar1.ln.intg=function(x=2, a_inv=0.5, mu=0.5, sig=2, y=1:3, u=c(0,1,1), v=c(0,2,2), s2=c(3,2,2))
{ #see ar1.intg; mu, sig : log normal parameters

  Pr=a_inv/(x+a_inv)
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem=c(tem,ar1.fun(y, u, v, s2, pr))
    }                                      
  res=tem*dlnorm(x, meanlog=mu, sdlog=sig) 
  return(res)
}

                                        # missing data dealt by the approximation
##ar1.fun oc to prob(Y=y|g)
ar1.fun=function(y=c(0,1,3), #counts
  U=c(0,1,1), # d*r[i,j-1]; r[i,j-1]=u[i,j]/alpha 
  V=c(0,2,2), # (1-d)*r[i,j-1]
  sz2=c(3,2,2), # r[i,j]- d*r[i, j-1]), 
  pr=0.5      # pr = prob =(1/a)/(g+1/a)
  )
{
  ## Compute: P(Y_i1=y_i1|G=g) P(Y_i1=y_i1|y_i2,G=g)...P(Y_ini=y_ini|y_i(ni-1),G=g)
  ## where P(Y_i1=y_i1|g_i) = NB(r_i1,pi)
  ##   and P(Y_ij=y_ij|y_i(j-1),g_i) = sum_{k=0}^{min{y_ij,y_ij-1} } P(Z=k|y_i(j-1),g_i)P(eps=y_ij-k|y_i(j-1),g_i)
  n=length(y) ## the number of repeated measure
  
  ## the first repeated measure P(Y_i1=y_i1|g_i)
  pb=dnbinom(y[1], size = sz2[1], prob=pr)
  if (n==1) return(pb)

  for ( i in 2:n)
    { ## P(Y_ij=y_ij|y_i(j-1),g_i)
      k = 0:min(y[i], y[i-1])
                                        #betabinomial dbb(x, N, u, v)
      pp=dbb(x=k, N=y[i-1], u=U[i], v=V[i])
      pp=pp*dnbinom(y[i]-k, size = sz2[i], prob=pr)
      p1=sum(pp)
      pb=pb*p1 
    }
  return(pb)
}
##============= end MLE ============

##============= semiparametric procedure (SP) ================
                                        #gwi = gi *wi + 1-wi 

mle.ar1.non3 <- function(
                         formula,     ## an object of class "formula"
                         ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
                         data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
                         ## containing the variables in the model.
                         ID,          ## a vector of length n*ni containing patient IDs of observations in data
                         Vcode,
                         p.ini=NULL, # log(a), log(var(G)), logit(d), b0, b1,... if its NULL then initial values are set to 0.1,NULL,0.1,log(mean(y))+0.1,0.1,0.1,...
                         IPRT=TRUE,
                         ##FixN=-1, #see mle.ar1.fun
                         deps=1.e-3,      # stop iteration when max(bhat.new - bhat.old) < deps
                         maxit=100
                         )
{
  dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat)
  
  ## DT = (id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames)

  ## DT$dif contains Vcode[j+1]-Vcode[j] except the initial scan of each patient (which is zero).
  DT$dif <- c(0,diff(Vcode))
  DT$dif[DT$ind] <- 0
  DT$dif <- DT$dif[1:DT$totN] #dif = scan lag 
  
  Dss <- getDjj(DT, Vcode)  # j-j' for all posible j > j'
  
  if (is.null(p.ini))
    {
      p.ini=rep(0, 4+DT$cn)
      p.ini[4]=mean(DT$y)   
    }
  p.new <- p.ini  #initial values for the first iteration
  W <- NULL       #weights for GLS

  dth <- 1  #dth = max(bhat.new-bhat.old)
  ahat <- 0 #initial value for alpha
  counter <- 1
  repeat
    {
      p.ini <- p.new #reset the initial values at each iteration
      b.ini <- p.ini[-(1:3)] ## b0,b1,...
      
      ## 
      ## Main Step 1 : estimate b using Generalized Least Square
      ## 
      if (DT$cn>0) #with covariates
        {
          temB <- optim(b.ini+0.1,## initial values of coefficients b0,b1,...
                        estb, ## weighted residuals function
                        hessian=TRUE,
                        dat=DT,## output of getDT
                        Wh=W,## list of N elements, containing ni by ni weight matrix for each patient
                        PRT=FALSE)
          bhat <- temB$par 
        }else{  #without covariates
          tem <- optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
          bhat <- tem$min
        }

      ## Stopping criterion
      ## check if max(bhat.new - bhat.old) < deps
      dth=max(abs(bhat-b.ini))
      ## if (IPRT) print(dth)
      if (dth < deps || counter > maxit) break

      ## 
      ## Main Step 2: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)
      ## 
      ## Substep 1: the estimation of mu_i+
      ##         hat.mu_i+ = exp( X_i^T bhat)
                                        # compute exp(b0+b1*x1+ ...)
      mu_ij <- rep(exp(bhat[1]), DT$totN)
      if (DT$cn>0) {
        tem <- exp(DT$x%*%bhat[-1]) 
        mu_ij <- tem*mu_ij ## vector of length Ntot containing hat.mu_ij
      } 
      ## mu_ip = u.i+ 
      mu_ip <- as.vector(tapply(mu_ij, DT$id, sum)) ## mu_i = sum_j mu_ij
      gh0 <- DT$ys/mu_ip ## vector of length n containing y_i+/mu_i+
      
      ## Substep 2: the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
      ## weights for gi
      ## compute wi =sqrt(var(G)/Var(Y+/u+))
      if (ahat == 0) #wi =1 at the first iteration
        {
          wi <- rep(1, DT$np)
        }else{
          gsig <- by(data=cbind(mu_ij, Vcode), INDICES=DT$id, FUN=getWhs.ar1,
                     th=that, a=ahat, dt=dht, act="sum")/mu_ip^2 
          wi <- sqrt(that/gsig)  #that =var(G), gsig=var(Y+/u+) 
        }
      
      gh <- wi*gh0+1-wi

      ## normalized gh so that the average is 1
      gh <- gh/mean(gh)
      gh <- as.vector(gh)
      
      if (IPRT) {
        cat("\n\n iteration",counter)
        cat("\n The distribution of hat.g \n")
        print(summary(gh))
      }

                                        # frequency table for gh
      tem <- sort(round(gh,6))
      gh1 <- unique(tem)
      ghf <- table(tem)
      gtb <- cbind(ghf, gh1)
      rownames(gtb) <- NULL
      gtb1 <- gtb
      gtb1[,1] <- gtb[,1]/DT$np #covert freq to prop
      ##estimate alpha and delta
      ## 
      ## Main Step 3: estimate a and d using profile likelihood
      ## 
      ad.ini=p.ini[c(1,3)]

      tt <- optim(ad.ini+0.1, #c(log(alpha), logit(delta))
                  ar1.ad.lk,
                  hessian=TRUE,
                  dat=DT,
                  Iprt=FALSE,
                  gTB=gtb1,#frequency table for hat(gi): with proportions
                  uhat=mu_ij #hat(u.ij)
                  ##,FixN=FixN
                  )
      
      ahat=exp(tt$par[1])
      dht=ilgt(tt$par[2])

      ## 
      ## Main Step 4: moment estimate of var(G)
      ##
      that=estSg3.ar1(DT, mu_ip, mu_ij, ahat, dht, Dss$dis) 
      
                                        #weights for GLS
      W <- by(data=cbind(mu_ij, Vcode), INDICES=DT$id,
              FUN=getWhs.ar1, th=that, a=ahat, dt=dht, act="inv") 

                                        #update the estimate
      p.new=c(tt$par[1], log(that), tt$par[2], bhat)

      if (IPRT){
        cat("\n log_a, log_var(G), lgt_d, log_u0", DT$xnames)
        cat("\n", p.new)
      }
      counter <- counter + 1
    }

  if (maxit==counter) warning("The maximum number of iterations occur!")
  vcm=NULL
  if (DT$cn>0) vcm=solve(temB$hessian) #vcm=vcm.fun(he=temB$hessian, corr=F)
  p.new <- matrix(p.new,ncol=1)
  rownames(p.new)=c( "log_a", "log_var(G)", "lgt_d", "(Intercept)", DT$xnames)
  if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c("(Intercept)", DT$xnames)
  re <- list(opt=tt, V=vcm, est=p.new, gi=as.vector(gh), mod="NoN",##idat = data.frame(dat),
             Vcode=Vcode,cor="ar1",formula=formula)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}



## getWhs.ar1=function(ud,  #cbind(u, sn)
##                     th=exp(1.3), #Var(G)
##                     a=exp(-0.5), dt=0.5, 
##                     act="inv" #  transformation option: "none", "inverse", "sum"
##                    )
## { u =ud[,1]
##   sn=ud[,2]

##   w1=u%*%t(u)*th # u.i*u.j*Var(G)

##   n=length(u)
##   Um=diag(u)%*%matrix(1,n,n)
##   mU=pmin(Um, t(Um)) #pairwise min (u.j, u.j')

##   Dm=diag(sn)%*%matrix(1,n,n)
##   mD=abs(Dm-t(Dm))  

##   Md=dt^mD #Mdis

##   UD=mU*Md #min(uj, uj')*d^(distance)
##   w2=(1+(th+1)*a)*UD   #[d^abs(j-j')]*u.l*(a*Var(G)+a+1); l=min(j,j')
  
##   W=w1+w2
##   if (act=="inv") W=solve(W)
##   if (act=="sum") W=sum(W)
##   return(W)
## }


getWhs.ar1 <- function(ud,  #cbind(mu_ij, Vcode)
                       th=exp(1.3), #Var(G)
                       a=exp(-0.5), dt=0.5, 
                       act="inv" #  transformation option: "none", "inverse", "sum"
                       )
{ mu_ij  <- ud[,1]
  Vcode <- ud[,2]

  w1 <- mu_ij%*%t(mu_ij)*th # u.i*u.j*Var(G)

  n <- length(mu_ij)
  Um <- diag(x=mu_ij,nrow=n)%*%matrix(1,n,n)
  mU <- pmin(Um, t(Um)) #pairwise min (u.j, u.j')

  Dm <- diag(x=Vcode,nrow=n)%*%matrix(1,n,n)
  mD <- abs(Dm-t(Dm))  

  Md=dt^mD #Mdis

  UD=mU*Md #min(uj, uj')*d^(distance)
  w2=(1+(th+1)*a)*UD   #[d^abs(j-j')]*u.l*(a*Var(G)+a+1); l=min(j,j')
  
  W=w1+w2
  if (act=="inv") W=solve(W)
  if (act=="sum") W=sum(W)
  return(W)
}



  ##                                       #example
##                                         #jk=ar1.exdt[ar1.exdt$Vcode<4,]
##                                         #mle.ar1.non3(dat=jk, p.ini=c(0,0,0, 0.5), FixN=5)

##                                         #no weights for gi
##                                         #mle.ar1.non(dat=jk[1:500,], p.ini=ar1.jk1fit$est[1:4,1])
##                                         #deps control for iteration stop when sum(abs(b[i]-b[i-1]))<deps
##                                         #see mle.ar1.non3
## mle.ar1.non <- function(dat, p.ini=NULL, IPRT=TRUE, ##FixN=-1,
##                         deps=1.e-3)
## {
##   DT <- getDT(dat[,-2])  
##                                         #DT=(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))
##   DT$dif <- c(0,diff(dat[,2]))
##   DT$dif[DT$ind] <- 0
##   DT$dif <- DT$dif[1:DT$totN]

##   Dss <- getDjj(DT, dat[,2]) 
  
##   if (is.null(p.ini))
##     { p.ini=rep(0, 4+DT$cn) 
##       p.ini[4]=log(mean(DT$y))
##     }
##   p.new=p.ini
##   W=NULL

##   dth=1 #max(p.new-p.ini)
##                                         #while(dth > 1.e-3)
##   repeat
##     { p.ini=p.new
##       b.ini=p.ini[-(1:3)]
      
##       if (DT$cn>0) 
##         { temB=optim(b.ini+0.1, estb, hessian=T, dat=DT, Wh=W, PRT=F)
##           bhat=temB$par }
##       else
##         { tem=optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
##           bhat=tem$min }

##       if (IPRT) print(bhat)
##       dth=max(abs(bhat-b.ini))
##       if (IPRT) print(dth)
##                                         #if (dth < 1.e-3) break
##       if (dth < deps) break

##       uh=rep(exp(bhat[1]), DT$totN)
##       if (DT$cn>0) { tem=exp(DT$x%*%bhat[-1])
##                      uh=tem*uh }

##       uhat=as.vector(tapply(uh, DT$id, sum))
##       gh0=DT$ys/uhat
      
##       gh_1=min(gh0[gh0>0])/2
##       ll=sum(gh0==0)
      
##       gh=gh0
##                                         #gh[gh0==0]=seq(0,gh_1, length=ll)
##       ll1=ll%%10
##       ll2=ll%/%10
##       ll3=rep(1:0, c(ll1, 10-ll1))+ll2
##       gh[gh0==0] = rep(seq(0, gh_1, length=10), ll3)


##       gh=gh/mean(gh)
##       gh=as.vector(gh)
      
##       tem=sort(round(gh,6))
##       gh1=unique(tem)
##       ghf=table(tem)
##       gtb=cbind(ghf, gh1)
##       rownames(gtb)=NULL

##       ad.ini=p.ini[c(1,3)]
##       tt=optim(ad.ini+0.1, ar1.ad.lk, hessian=T, dat=DT, Iprt=IPRT, gTB=gtb, uhat=uh##, FixN=FixN
##         )
      
##       ahat=exp(tt$par[1])
##       dht=ilgt(tt$par[2])

##       that=estSg3.ar1(DT, uhat, uh, ahat, dht, Dss$dis)
##       W=by(cbind(uh, dat[,2]), DT$id, getWhs.ar1, th=that, a=ahat, dt=dht, act="inv") 
##       p.new=c(tt$par[1], log(that), tt$par[2], bhat)

##                                         #dth=max(abs(p.new-p.ini))
##       if (IPRT) print(p.new)
##     }
  
##                                         #nlk=tt$value+sum(lgamma(DT$y+1))
##   vcm=NULL
##   if (DT$cn>0) vcm=solve(temB$hessian) #vcm.fun(he=temB$hessian, corr=F)
##                                         #get V(b)
##                                         #var(b)=(X'D(V*)^{-1}DX)^{-1}
##                                         #(V*)=th+diag((1+a(th+1))/uij)
##                                         #p.est=c(log(uhat), tt$min, log(that))
##   names(p.new)=c( "log_a", "log_th", "lgt_d", "log_u0", DT$xnames)
##                                         #if (IPRT) print(p.new)
##   return(list(opt=tt, vcm=vcm, est=p.new, gi=as.vector(gh), mod="NoN"))
## }



##subroutines for mle.ar1.non1 and mle.ar1.non3
                                        #likelihood function for estimating alpha and delta
ar1.ad.lk <- function(para=c(-0.5, -1), #c(log(alpha), logit(delta)) 
                      dat,  #see ar1.lk
                      Iprt=TRUE, 
                      gTB, #frequency table for hat(gi): with proportions
                      uhat #hat(u.ij)
                      ##,FixN=9  #see ar1.lk
                      )
{
  ## Li(a,d;bhat,ghat,ys)
  ## = 1/N* sum_{l=1}^N Pr(Y_i,1=y_i,1;ghat,bhat,a)*prod_{j=2}^ni Pr(Y_ij=y_ij;y_i(j-1),ghat,bhat,a) 
  if (Iprt) cat(para)
  ainv=exp(-para[1]) #ainv=1/a  
  dt=ilgt(para[2])   #delta
  
  th2=uhat*ainv  #r.ij
  Dl=dt^dat$dif  

  szm=c(0, th2[-dat$totN])
  U=Dl*szm
  V=szm-U
  sz2=th2-U
  tem=dat$ind[1:dat$np]
  sz2[tem]=th2[tem]
  if (any(sz2<=0)) return(1.e15)

  nllk=0
  lk0=0 #likelihood for a patient with all 0 counts
  
  Pr=ainv/(ainv+gTB[,2]) ## in this manuscript p=p
                                        #fq=gTB[,1]/dat$np      #proportions
  
  for (i in 1:dat$np)
    {

      ##if (lk0>0&dat$ni[i]==FixN&dat$ys[i]==0) lki=lk0
      ##else{
      ## Li(a,d;bhat,ghat,ys)
      ## = 1/N* sum_{l=1}^N Pr(Y_i,1=y_i,1;ghat,bhat,a)*prod_{j=2}^ni Pr(Y_ij=y_ij;y_i(j-1),ghat,bhat,a)
      ll=dat$ind[i]:(dat$ind[i+1]-1)
      tem=ar1.non(Pr=Pr, y=dat$y[ll], u=U[ll], v=V[ll], s2=sz2[ll])
      lki=sum(tem*gTB[,1])
      ##if (dat$ys[i]==0&dat$ni[i]==FixN) lk0=lki
      ##}
      nllk = nllk - log(lki)
    }
  if (Iprt) cat(" nllk=", nllk, "\n")
  return(nllk)
}

                                        #ar1.non=function(Pr=Pr, y=1:3, sz=rep(3,3), dt=0.3)
ar1.non=function(Pr, y=1:3, u=c(0,1,1), v=c(0, 2,2), s2=c(3,2,2))
{ #compute ar1.fun for a vector of Pr
                                        #also see ar1.intg
  tem=NULL
  for ( pr in  Pr)
    { #print(pr)
      tem=c(tem,ar1.fun(y, u, v, s2, pr)) 
    }
  return(tem)
}


                                        #moment estimate of var(G)
                                        #estSg3.ar1(DT, uhat, uh, ahat, dht, Dss$dis)
estSg3.ar1=function(dat, #see ar1.lk
  ui, uij,     #u.i+ and u.ij
  ahat, dhat,  #alpha hat and delta hat
  Dis          # j - j'
  ) 
{ 
  ssy=sum(dat$ys^2)-sum(dat$y^2)
  ssu=sum(ui^2)-sum(uij^2)
  
  DU=2*sum(dhat^Dis[,1]*uij[Dis[,2]]) 
  that=(ssy-DU)/(ssu+DU*ahat)
                                        #print(ssy/ssu)
  return(max(that-1, 0))
}

                                        # compute j-j' for all posible j > j'
getDjj <- function(DT,Vcode)
{
  kk <- Nj <- dj <- NULL
  for (i in 1:DT$np)
     {
       if (DT$ni[i] > 1)
         {
           ll <- DT$ind[i]:(DT$ind[i+1]-1)
           iVcode <- Vcode[ll]
           nj <- djj <- NULL
           id <- DT$id[DT$ind[i]]
           
           Ni <- DT$ni[i]-1
           for (j in 1:Ni )
             { 
               for (jj in (j+1):DT$ni[i])
                 djj <- c(djj, iVcode[jj]-iVcode[j])
             }
           nj <- rep(1:Ni,Ni:1)
           kk <- c(kk,rep(id, length(nj)))
           Nj <- c(Nj, nj+DT$ind[i])
           dj <- c(dj, djj)
         }
     }
   return(list(id=kk, dis=cbind(dj, Nj-1)))
                                        # dis=([i]j'-[i]j, index for u.ij)
 }

##route for weights
## act= inv => W = Var(Y)(^-1)
## act= sum => W= sum(Var(Y))
## Var(Y)=Cov(Yj, Yj')=u.j*u.j'*Var(G)+[d^abs(j-j')]*u.l*(a*Var(G)+a+1); 
## where l=min(i,j)

                                        #example ud.ex=cbind(rep(exp(.48),5), 1:5)
                                        #getWhs.ar1(ud=ud.ex, th=exp(1.46), a=exp(-0.825), dt=ilgt(-0.26), act="sum")
##==========  end SP =============

##==========  Index functions  ============
                                        # pr(q(Ynew) >= q(ynew) | Ypre=ypre)
                                        # input: data from one patient

jCP.ar1=function(tpar, ## log(a),log(theta),log(delta),b0,...
                                        # parameters (obj is an output from mle.ar1.non3 or mle.ar1.fun) 
  ypre,## vector of length # pre, containing CEL
  ynew,  ## vector of length # new, containing CEL
  y2m=NULL, 
  XM, # matrix of covariates
  stp,  
  mod="G", #options for RE dist'n: G=Gamma, N=lognormal, NoN=nonparametric
  LG=FALSE,    # if T, return logit(P)
  MC=FALSE, N=40000, #Monte carlo integration 
  qfun="sum", #q function
  oth=NULL    #if mod=="NoN", oth=obj$gi
  )
{  
  a = exp(tpar[1])
  th = exp(tpar[2])
  dt = ilgt(tpar[3])
  
  sn=length(ypre)+length(ynew) ## the number of repeated measures 
  
  u0=exp(tpar[4]) ## beta0
  u=rep(u0, sn)
  ## with covariates
  if (length(tpar)>4) u=u*exp(XM%*%tpar[-(1:4)]) 
  
  if (MC)
    tem=MCCP.ar1(ypre=ypre, ynew=ynew, stp=stp, u=u, th=th, a=a, dt=dt, mod=mod, Ns=N, gh=oth, qfun=qfun)
  else 
    tem=CP1.ar1(ypre=ypre, ynew=ynew, y2m=y2m, stp=stp, u=u, th=th, a=a, dt=dt, mod=mod, gh=oth, qfun=qfun)
  if (LG) tem=lgt(tem)
  return(tem)
}

##subroutines for jCP.ar1
                                        #CP1.ar1 (copied Psum.jun25.R)
                                        #parameters on the original scales

                                        #CP1.ar1(ynew=c(1,0,1), mod="NoN", gh=obj$gi, qfun="max")
CP1.ar1=function(ypre=c(0,1),
  ynew=c(1,0,1),
  y2m=NULL, 
  stp=c(0, 1, 1, 1, 1), 
  u=rep(1.5,5),        
  ## sn by # covariate matrix containing:
  th=3, a=0.5, dt=1/3, #parameters on the original scales
  mod="G",  #G=gamma, NoN=nonparam
  gh, ## NULL for parametric model
  qfun="sum")
  {
    if (all(is.na(ynew))) return (NA)
    
    Qf=match.fun(qfun)
    newQ=Qf(ynew, na.rm=T)
    if (newQ==0) return(1)

    ain=1/a
    if (mod=="G") shp=1/th
    if (mod=="NoN") 
      { tem=sort(round(gh,6))
        gh1=unique(tem)
        gp=table(tem)/length(tem)
                                        #Pr=ain/(ain+gh1)
        Pr=1/(gh1*a+1)
      }
    
    n0=length(ypre)
    n1=length(ynew)
    l1=n0+n1  

    y=c(ypre, ynew) 
                           

    sz=u/a
    DT=dt^stp
    
                                        #parameters for ypre
    szm=c(0, sz) 
    szm=szm[-length(szm)]
    u=DT*szm
    v=szm-u
    sz2=sz-u                           #possible combinations for ynew under q()
                                        #y2m=getY2.1(newQ, l1-l0, qfun)
    if (is.null(y2m)) y2m=getY2.1(newQ, l1-n0, qfun) 
    ## l1-n0 = n1 # the number of new-scans
    ## newQ = q(Y_new) 
    if (mod=="G")
      { #Pr(Ypre=Ypre)
        tem=integrate(ar1.intg, lower=0, upper=Inf,  
          a_inv=ain, sh=shp, sc=th, 
          y=ypre, u=u[1:n0], v=v[1:n0], s2=sz2[1:n0])
        bot=tem$v
                                        #print(bot)
        tem=integrate(ar1.2.mmg, lower=0, upper=Inf, rel.tol=0.1, 
          a_inv=ain, sh=shp, sc=th, 
          y1=ypre, y2m=y2m, u=u, v=v, sz2=sz2)
        top=tem$v
                                        #print(top)
      }
    if (mod=="N") #added May 18, 2012
      {
        ## YK: June 6 {
        s.ln=log(th+1)    
        u.ln=-s.ln/2     
        s.ln=sqrt(s.ln)
        ## }
        tem=integrate(ar1.ln.intg, lower=0, upper=Inf,  
          a_inv=ain, mu=u.ln, sig=s.ln, 
          y=ypre, u=u[1:n0], v=v[1:n0], s2=sz2[1:n0]) 
        bot=tem$v
        tem=integrate(ar1.ln2.mmg, lower=0, upper=Inf, rel.tol=1.e-3, 
          a_inv=ain, mu=u.ln, sig=s.ln, 
          y1=ypre, y2m=y2m, u=u, v=v, sz2=sz2)
        top=tem$v
      }

    if (mod=="NoN")
      { p1.non=p2.non=NULL
        for ( pr in Pr)
          { #p1=ar1.fun(y1, u[1:l0], v[1:l0], sz2[1:l0], pr)
            p1=ar1.fun(ypre, u[1:n0], v[1:n0], sz2[1:n0], pr)

            p1.non=c(p1.non, p1)
            
            p2=0
            for (j in 1:nrow(y2m))
              { #yy=c(y1[l0], y2m[j,]) 
                                        #p2=p2+ar22.fun(yy, u[l0:l1], v[l0:l1], sz2[l0:l1], pr) 
                yy=c(ypre[n0], y2m[j,]) 
                p2=p2+ar22.fun(yy, u[n0:l1], v[n0:l1], sz2[n0:l1], pr) 
              }
            p2.non=c(p2.non, p2*p1) 
          }
        bot=sum(p1.non*gp)
        top=sum(p2.non*gp)
      }

                                        #print(c(bot, top))
    res=1-top/bot
    return(res)
  }


                                        # Pr( Ypre=y1, Ynew[1:k-1] = y2m[i,1:(k-1)], Ynew[k]<=y2m[i,k] | G=x) 
                                        # with G~gamma for all row of y2m
ar1.2.mmg=function(x=0.5, a_inv=2, sh=1/3, sc=3, y1=c(0,1), 
  y2m=getY2.1(), #a matrix, each row is a set of value for Ynew
  u=c(0,1,1,1,1), v=c(0,2,2,2,2), sz2=c(3,2,2,2,2))
{  #print(x)
                                        #print(length(x))
  
  m=nrow(y2m)
  Pr=a_inv/(a_inv+x)
  l0=length(y1)
  l1=ncol(y2m)+l0     

  tem=NULL
  for ( pr in Pr)
    { #Pr(Ypre=ypre)
      
      p1=0
      for (j in 1:m)
        {  
          y=c(y1[l0], y2m[j,]) 
          p1=p1+ar22.fun(y=y, U=u[l0:l1], V=v[l0:l1], sz2=sz2[l0:l1], pr=pr) 
        }
      p1=p1*ar1.fun(y=y1, U=u[1:l0], V=v[1:l0], sz2=sz2[1:l0], pr=pr) 
      tem=c(tem, p1)
    }
  tem=tem*dgamma(x, shape=sh, scale=sc)
  return(tem)
}

ar1.ln2.mmg=function(x=0.5, a_inv=2, mu=1/3, sig=3, y1=c(0,1), 
  y2m=getY2.1(), #a matrix, each row is a set of value for Ynew
  u=c(0,1,1,1,1), v=c(0,2,2,2,2), sz2=c(3,2,2,2,2))
{  
  m=nrow(y2m)
  l0=length(y1)
  l1=l0+ncol(y2m)
  Pr=a_inv/(a_inv+x)
  
  tem=NULL
  for ( pr in Pr)
    { #Pr(Ypre=ypre)
      
      p1=0
      for (j in 1:m)
        {   
          y=c(y1[l0], y2m[j,])
          p1=p1+ar22.fun(y=y, U=u[l0:l1], V=v[l0:l1], sz2=sz2[l0:l1], pr=pr) 

        }
      p1=p1*ar1.fun(y=y1, U=u[1:l0], V=v[1:l0], sz2=sz2[1:l0], pr=pr)
      tem=c(tem, p1)
                                        #print(pr)
    }
  tem=tem*dlnorm(x, meanlog=mu, sdlog=sig)
  return(tem)
}

                                        #Pr(Y2=y2, Y3=y3, ..., Yn<= yn|Y1=y1)
ar22.fun=function(y=c(0,1,3), U=c(0,1,1), V=c(0,2,2), sz2=c(3,2,2), pr=0.5)
{  n=length(y)   

   pb=1
   for ( i in 2:n)
     { k = 0:min(y[i], y[i-1])
                                        #betabinomial dbb(x, N, u, v)
                                        #source("/home/yinshan/Mydesk/Mylib/Rlib/myfun/beta-binomial.R")
       pp=dbb(x=k, N=y[i-1], u=U[i], v=V[i])
       if (i<n) { pp=pp*dnbinom(y[i]-k, size = sz2[i], prob=pr)}
       else  { pp=pp*pnbinom(y[i]-k, size = sz2[i], prob=pr) }
       p1=sum(pp)
       pb=pb*p1 
     }
   return(pb)
 }


                                        #MCCP.ar1: monto carlo for conditional probability
                                        #MCCP.ar1=function(ypre, ynew, u, th, a, dt, mod="G", Ns=1000, gh=gi, qfun=sum)
MCCP.ar1=function(ypre, ynew, stp, u, th, a, dt, mod="G", Ns=1000, gh, qfun="sum") 
{ if (all(is.na(ynew))) return (NA)

  Qfun=match.fun(qfun)
  newQ=Qfun(ynew, na.rm=T)
  if (newQ==0) return(1)

  if (newQ==1) 
    {  res=CP1.ar1(ypre=ypre, ynew=ynew, stp=stp, u=u, th=th, a=a, dt=dt, mod=mod, gh=gh,qfun=qfun)
       return(res)
     }

  ain=1/a
  if (mod=="G") shp=1/th

  if (mod=="N") 
    { 
      s.ln=log(th+1)   
      u.ln=-s.ln/2     
      s.ln=sqrt(s.ln)  
    }

  if (mod=="NoN") 
    { tem=sort(round(gh,6))
      gh1=unique(tem)
      gp=table(tem)/length(tem)
      Pr=1/(1+a*gh1)
    }
  
  n0=length(ypre)
  n1=length(ynew)
  l1=n0+n1  

  sz=u/a

  dt=dt^stp 
  szm=c(0, sz)
  szm=szm[-length(szm)]
  u0=dt*szm 
  v0=szm-u0
                                        #sz2=sz0[ms0]-u0  #YZ remove (May 24, 2012)
  sz2=sz-u0   
  Opre=list(y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0])
  Onew=list(y0=ypre[n0], u=u0[-(1:n0)], v=v0[-(1:n0)], s2=sz2[-(1:n0)], n1=n1)
                                        
  
  y=c(ypre, ynew) 
  if (mod=="G")
    { #Pr(Ypre=Ypre)
      tem=integrate(ar1.intg, lower=0, upper=Inf,  
        a_inv=ain, sh=shp, sc=th, 
                                        #y=ypre[ms0], u=u0, v=v0, s2=sz2) #YZ old
        y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0]) #YZ new (May 24, 2012)
      bot=tem$v

      tem=integrate(ar1.mmg, lower=0, upper=Inf, rel.tol=1.e-2, 
        a_inv=ain, sh=shp, sc=th, 
        O1=Opre, O2=Onew, tot=newQ, N=Ns, qfun=qfun) #YZ new (May 24, 2012)
                                      
      top=tem$v
    }
                                       
  if (mod=="N")
    { #Pr(Ypre=Ypre)
      tem=integrate(ar1.ln.intg, lower=0, upper=Inf,  
        a_inv=ain, mu=u.ln, sig=s.ln, 
        y=ypre, u=u0[1:n0], v=v0[1:n0], s2=sz2[1:n0])
      bot=tem$v

                                        #Pr(Ypre=ypre, Ynew+ < newSum) 
      tem=integrate(ar1.mm.ln, lower=0, upper=Inf, rel.tol=1.e-2, 
        a_inv=ain, mu=u.ln, sig=s.ln, 
        O1=Opre, O2=Onew, tot=newQ, N=Ns, qfun=qfun) 
      
      top=tem$v
    }
                            

  if (mod=="NoN")
    { p1.non=p2.non=NULL
      for ( pr in Pr)
        { 
          p1=ar1.fun(ypre, u0, v0, sz2, pr)  
          p2=p1*mmS.fun(Obj=Onew, pr=pr, nSum=newQ, N=Ns, qfun=qfun)

          p1.non=c(p1.non, p1)
          p2.non=c(p2.non, p2)
        }
      bot=sum(p1.non*gp)
      top=sum(p2.non*gp)
    }

                                        #print(c(bot, top))
  res=1-top/bot
  return(res)
}

##Pr(q(Ynew) < nSum |y0; pr) using MC
                                       
mmS.fun=function(Obj=list(y0=0, u=c(1,1,1), v=c(2,2,2), s2=c(2,2,2), n1=3),
  pr=0.5, nSum=2, N=1000, qfun="sum")
{ Qfun=match.fun(qfun)

                                        #k=length(sz1)
                                        #u=sz1[-k]*dt
                                        #v=sz1[-k]-u
                                        #s2=sz1[-1]-u

  u=Obj$u
  v=Obj$v 
  s2=Obj$s2

  y0=Obj$y0
  YY=NULL
  for (i in 1:Obj$n1)
    {  
      z1=rbb(N, y0, u[i], v[i])
      z2=rnbinom(N, size=s2[i], prob=pr)
      y1=z1+z2
      YY=cbind(YY, y1)
      y0=y1
    }
  if (Obj$n1>1) 
    { tot=apply(YY, 1, Qfun) 
    }
  else tot=YY

  return(mean(tot<nSum))
}


ar1.mmg=function(x, a_inv, sh, sc, O1, O2, tot, N, qfun=sum)  
{  #print(x)
                                        #print(length(x))
  Pr=a_inv/(a_inv+x)
  tem=NULL
  for ( pr in Pr)
    { #Pr(Ypre=ypre)
      p1=ar1.fun(y=O1$y, U=O1$u, V=O1$v, sz2=O1$s2, pr=pr)   
      p2=mmS.fun(Obj=O2, pr=pr, nSum=tot, N=N, qfun=qfun) 
                                        #print(c(pr, p1, p2))
      tem=c(tem,p1*p2)
    }
  tem=tem*dgamma(x, shape=sh, scale=sc)
  return(tem)
}

                                       
ar1.mm.ln <- function(x, a_inv, sig, mu, O1, O2, tot, N, qfun=sum)
{  Pr=a_inv/(a_inv+x)
   tem=NULL
   for ( pr in Pr)
     {
       p1=ar1.fun(y=O1$y, U=O1$u, V=O1$v, sz2=O1$s2, pr=pr)   
       p2=mmS.fun(Obj=O2, pr=pr, nSum=tot, N=N, qfun=qfun)
       tem=c(tem,p1*p2)
     }
   tem=tem*dlnorm(x, meanlog=mu, sdlog=sig)
   return(tem)
 }



dbb <- function(x, N, u, v) {
  beta(x+u, N-x+v)/beta(u,v)*choose(N,x)
}

## pbb <- function(q, N, u, v) {
##   sapply(q, function(xx) sum(dbb(0:xx, N, u, v)))
## }

## qbb <- function(p, N, u, v) {
##   pp <- cumsum(dbb(0:N, N, u, v))
##   sapply(p, function(x) sum(pp < x))
## }

rbb <- function(n, N, u, v) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}

                                        #logit and inverse logit
lgt <- function(p)
{
  log(p/(1-p))
}

ilgt <- function(x)
{
  tem=exp(x)
 res=tem/(1+tem)
 return(res)
}







CP.ar1.se=function(
  tpar,
  ypre,
  ynew,
  y2m=NULL,
  XM, 
  ## see jCP.ar1
  stp, 
  dist="G",   # dist'n of REs G=Gamma, N =lognormal
  V,  
  mc=FALSE,       # if true use MC
  qfun="sum")
{  if (all(is.na(ynew))) return(c(NA, NA))
   Qf=match.fun(qfun)
   newQ=Qf(ynew, na.rm=T)
   if (newQ==0) return(c(1,0))
   
   p=jCP.ar1(tpar=tpar,## log(a),log(theta),log(delta),b0,...
     ypre=ypre,## vector of length # pre, containing CEL
     ynew=ynew,## vector of length # new, containing CEL
     y2m=y2m, stp=stp,
     XM=XM, mod=dist, MC=mc, qfun=qfun,LG=FALSE)
   
   if (mc) mth="simple" else mth="Richardson"
                                       
   jac=jacobian(func=jCP.ar1, x=tpar, method=mth,
     method.args=list(eps=0.01, d=0.01, r=2), ypre=ypre, ynew=ynew,
     y2m=y2m, XM=XM, stp=stp, mod=dist, LG=TRUE, MC=mc, qfun=qfun)  
   s2=jac%*%V%*%t(jac)
   s=sqrt(s2) #SE of logit(Phat)
   return(c(p,s)) ## (Phat, SE(logit(Phat)))
 }
