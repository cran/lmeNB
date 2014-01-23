index.batch <- function(
                        data,
                        labelnp,
                        ID,
                        Vcode,
                        olmeNB=NULL, #model fit,
                        subset=NULL, 
                        ## compute the index for a subset of patients when subset = a vector of patient ids.
                        ## by default (subset=NULL), compute the index for all patients.
                        ## The indecies of ID must agree to the indices of idat0$ID, returned from olmeNB 
                        qfun="sum", # sum or max
                        IPRT=TRUE, #print control
                        i.se=TRUE, 
                        iMC=FALSE #for ar1 only
                        ## If olmeNB  ==  NULL then this is required
                        ){
  if (nrow(data) !=length(labelnp))stop("the length of labelnp does not agree with the number of rows in data")
  if (nrow(data) !=length(ID))stop("the length of labelnp does not agree with the number of rows in data")
  ## If the labelnp is TRUE/FALSE then the following code change them to 1/0
  ## If the labelnp is 1/0 then the following code keep them as 1/0
  labelnp[labelnp] <- 1; labelnp[labelnp] <- 0
  
  ftd <- formulaToDat(formula=olmeNB$formula,data=data,ID=ID,labelnp=labelnp)
  fulldata0 <- data.frame(ftd$dat)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  lnp <- ftd$labelnp
  if (olmeNB$cor != "ar1") Vcode <- rep(NA,length(fulldata0$ID))
  fulldat <- data.frame(ID=fulldata0$ID,Vcode=Vcode,CEL=fulldata0$CEL)
  ## If there is covariate, then add them in idat
  if (ncol( fulldata0 ) > 2) fulldat <- data.frame(fulldat,x=fulldata0[,3:ncol(fulldata0)])
  ## idat = dataframe(ID, Vcode, CEL, x.1, x.2, ...)

  if (is.null(subset)) id <- unique(fulldat$ID)  else id <- subset
  
  est <- as.data.frame(olmeNB$est)[,1] 
  n <- length(id)
  ## PP is output of this function
  PP <- matrix(NA, nrow=n, ncol=4)
  colnames(PP) <- c("p", "SE of logit(p_hat)", "low95", "up95")
  rownames(PP) <- unique(ftd$origID)

  l <- 0
  for ( i in id ) ## 1 to N (the total number of patietns)
    {
      l <- l+1
      ## extract n_i by ncol(fulldat) observations, corresponding to the patient i
      idat <- fulldat[fulldat$ID==i,,drop=FALSE] 
      
      if (nrow(idat)==0) next  #next if no data from this patient
      ilnp=lnp[fulldat$ID==i]
      
      if (! (1 %in% ilnp)) next  ## next if no new scans
      ##if (! (0 %in% ilnp)) next   ## next if no old scans !!!!!!CAUSUE problem for AR/semipara.

      xm <- NULL
      if ( ncol(idat)>3) xm <- as.matrix(idat[,-c(1:3)])
      ## the component of ilnpTF = TRUE if the repeated measure is pre scans.
      ilnpTF <- (ilnp==0)       

      if (olmeNB$cor=="ind")  ## independent model
        {
          ## total count on pre scans
          Y1 <- sum(idat$CEL[ilnpTF]) 
          ## if qfun = "sum" then Y2 = sum(idat$CEL[!ilnpTF]) else if qfun="max" then Y2= max(idat$CEL[!ilnpTF])
          Y2 <- eval(call(qfun, idat$CEL[!ilnpTF]))
                                      
          sn1 <- sum(ilnpTF)        ## number of pre scans
          sn2 <- sum(!ilnpTF)       ## number of new scans

          if (olmeNB$mod == "NoN"|i.se == FALSE) {
              ## If nonparametric method is used for the random effects
              tem <- jCP(tpar=est, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, 
                         XM=xm, dist=olmeNB$mod, type=qfun, oth=olmeNB$gtb)    
            }else{
              ## If distributional assumption was made for random effects
              tem <- CP.se(tpar=est, ##  log_a     log_th     log_u0
                           Y1=Y1, ## the sum of the response counts in pre
                           Y2=Y2, ## q(y_i,new)
                           sn1=sn1,## the number of pre scans
                           sn2=sn2,## the number of old scans
                           XM=xm, ## ni by # covariates matrix or NULL if no covariate
                           V=olmeNB$V,## covariance matrix of the maximum likelihood estimators: (log(a), log(theta), beta0, beta1, ...)
                           dist=olmeNB$mod,
                           pty=qfun)
            }  
        }else{          ##AR(1) model
                                     
          ypre <- idat$CEL[ilnpTF]
          ynew <- idat$CEL[!ilnpTF]
          stp <- c(0,diff(idat$Vcode))             
          Qf <- match.fun(qfun)
          newQ <- Qf(ynew, na.rm=T)
          ## possible combinations for ynew under q()
          y2m <- getY2.1(newQ, sum(!is.na(ynew)), qfun)
          ## the integration method is forced to be MC if the number of combinations that returns ynew under q() is large
          if (nrow(y2m) > 2000) 
            { iMC=TRUE
              y2m=NULL
            }
          
          if (olmeNB$mod =="NoN"|i.se==F) 
            {  
              tem <- jCP.ar1(tpar=est, ypre=ypre, ynew=ynew, 
                             y2m=y2m, 
                             stp=stp, 
                             XM=xm,   
                             mod=olmeNB$mod, MC=iMC, qfun=qfun, oth=olmeNB$gi)

            }else{            
              tem <- CP.ar1.se(
                               tpar=est, ## log(a),log(theta),log(delta),b0,...
                               ypre=ypre, ## vector of length # pre, containing CEL
                               ynew=ynew, ## vector of length # new, containing CEL
                               y2m=y2m, 
                               XM=xm,  
                               stp=stp,
                               dist=olmeNB$mod, V=olmeNB$V, mc=iMC, qfun=qfun)
            }
        }
      if (olmeNB$mod == "NoN" | i.se==FALSE) {  
        tem1= c(NA, NA, NA) 
      }else{ tem1=pp.ci(tem) }
      PP[l,]=c(tem, tem1)
      
      if (IPRT)
        {
          round.tem <-round(tem,3)
          if (olmeNB$mod != "NoN" & i.se==TRUE) 
            { 
              round.tem1 <- round(tem1,3) 
              cat("\n patient", rownames(PP)[i]) 
              cat(": estimated Pr (SE of logit(p_hat)) ",round.tem[1]," (",round.tem[2],")", sep="")
              cat( ": 95% CI [",round.tem1[1], "\ ", round.tem1[2], "]\n", sep="")
            }else{
              cat("\n patient",  rownames(PP)[i], "estimated Pr",round.tem)
            }
        }
    }  
  if (IPRT) cat("\n") 


  res <- list()
  res$condProbSummary <- PP
  res$para$CEL <- model.response(model.frame(formula=olmeNB$formula,data=data,na.action=function(z)z)) ## NA is kept
  res$para$labelnp <- labelnp ## original input which contains labelnp for corresponding NA
  res$para$ID <- ID ## original ID which contains ID for the corresponding NA
  res$para$qsum <- qfun
  class(res) <- "IndexBatch"
  return(res)
}

getY2.1=function(tot=2, ## q(Y_new)
  n=3, ## # new scans
  qfun="sum")
{
  if (tot==0){ return(matrix(rep(0,n),nrow=1 ))}
  if (n==1) 
    { ym=tot-1
      return(matrix(ym, ncol=n))
    }

  N=tot^(n-1)
  ii=0:(N-1)
  
  ym=NULL
  for (j in 1:(n-1))
    { t2=ii%%tot
      ii=ii%/%tot
      ym=cbind(ym, t2)
    }
  if (qfun=="sum") 
    { ysum=apply(ym, 1,sum)
      yn=tot-ysum-1
      ym=cbind(ym, yn)
      ym=ym[yn>=0,]
    }
  else { yn=rep(tot-1, nrow(ym)) 
         ym=cbind(ym, yn) }
  
  return(matrix(ym, ncol=n))
}



Psum1 <-
  function(Y1, Y2,     #Y1=sum(y.pre), Y2=sum(y.new)
           u1, u2, # u1=mean(Y1), u2=mean(Y2)
           a, Th,           #Var(Gi), under the gamma model, Th=scale
           dist="G",       #"G" = gamma, "N" = log-normal, "GN" = mix of gamma and normal#"U" = log-uniform #"NoN" = non-parametric
           othr=NULL,       # if dist= "GN",
           sn1
           ## othr=list(u.n=3, s.n=0.5, p.mx=0.05, sh.mx=NA)
           ## if dist="NoN", 
           ## othr = ghat frequency table of gi (olmeNB$gtb)
           ## for other dist options, othr = NULL 
           )
  { if (Y2==0) {return(1)}
    uVa=c(u1, u2)/a
    if (dist=="NoN")
      { 
        tem=Psum.non(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, gi=othr)
        return(tem)
      }
    else if (dist=="G") #gamma
      { 
        tem1 <- integrate(cum.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, y1=Y1, y2=Y2, u1=uVa[1], u2=uVa[2], abs.tol=1e-75,sn1=sn1)
        if (sn1 > 0) tem2 <- integrate(int.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, ysum=Y1, usum=uVa[1], abs.tol=1e-75)
        else{
          tem2 <- list();
          tem2$v <- 1
        }
      }  
    else if (dist=="N")
      {  t1 <- sqrt(log(Th+1)) #sd (sc)
         t2 <- -t1^2/2 #mean (sh)
         tem1 <- integrate(cum1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, y1=Y1, y2=Y2, u1=uVa[1], u2=uVa[2], abs.tol=1e-75,sn1=sn1)
         if (sn1 > 0) tem2 <- integrate(int1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, ysum=Y1, usum=uVa[1],abs.tol=1e-75)
         else{
           tem2 <- list();
           tem2$v <- 1
         }
       }
    else {
      stop("mod must be G, N or NoN")
    }

    val=tem1$v/tem2$v
    return(1-val)
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

Pmax1 <-
  function(Y1=0,                #sum(ypre)
           Y2=1,                #max(ynew)
           u1=3,                #Ex(sum(Ypre))
           u2=c(1.5, 1.5, 1.5), #Ex(Ynew); vector 
           a=0.5, 
           Th=3,                #var(G), no use if dist = "NoN" 
           dist="G",            # "G"=gamma, "N" = lognormal, "NoN = nonparametric
           othr=NULL,            # othr=gi (a vector) if  dist = "NoN", not used otherwise
           sn1
           )
{ if (Y2==0) {return(1)}
  uVa1=u1/a
  uVa2= u2/a
  
  if (dist=="NoN")
    { tem=Pmax.non(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, gi=othr)
      return(tem)
    }
  else if (dist=="G") 
    {
      tem1 <- integrate(max.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, y1=Y1, maxY2=Y2, u1=uVa1, u2=uVa2, abs.tol=1e-75,sn1=sn1)
      if (sn1 > 0) tem2 <- integrate(int.fun, lower=0, upper=Inf, a_inv=1/a, sh=1/Th, sc=Th, ysum=Y1, usum=uVa1, abs.tol=1e-75)
      else{
        tem2 <- list();tem2$v <- 1
      }
    }  
  else if (dist=="N")
    {  t1=sqrt(log(Th+1)) #sd (sc)
       t2=-t1^2/2 #mean (sh)
       tem1 <- integrate(max.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, y1=Y1, maxY2=Y2, u1=uVa1, u2=uVa2, abs.tol=1e-75,sn1=sn1)
       if (sn1 > 0) tem2 <- integrate(int1.ln, lower=0, upper=Inf, a_inv=1/a, sh=t2, sc=t1, ysum=Y1, usum=uVa1,abs.tol=1e-75)
       else{
         tem2 <- list();tem2$v <- 1
       }
     }
  else stop("mod must be G, N or NoN")
  val <- tem1$v/tem2$v
  return(1-val)
}


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



## ================ numerator of CPI ===========================

cum.fun <-
  function(
           ## Pr(Y_i,new+ >= y_i,new+, Y_i,pre+ = y_i,pre+, G_i=g) 
           ## = Pr(Y_i,new+ >= y_i,new+|Y_i,pre+ = y_i,pre+, G_i=g) Pr(Y_i,pre+ = y_i,pre+|G_i=g) Pr(G_i=g)
           ## = Pr(Y_i,new+ >= y_i,new+|G_i=g) Pr(Y_i,pre+ = y_i,pre+|G_i=g) Pr(G_i=g) ## conditional independence
           ## = (1-pnbinom(Y_i,new+;size=u2,prob=p))*dnbinom(Y_ipre+;size=u1,prob=p)*dgamma(gi;shape=sh,scale=1/sh)
           ## where u1 = sum_{j in old scan} exp(beta^T * Xij) and u2 = sum_{j in new scan} exp(beta^T * Xij)  
           x=2,          # value of the random effect
           a_inv=0.5,    # 1/a
           sh=0.5, sc=2, # shape and scale of the Gamma dist'n
           y1=2, y2=2,   # Y1=sum(y.pre), Y2=sum(y.new)
           u1=3, u2=3,    # u1 = r1 = mu1/a; u2= r2 =mu2/a
           sn1
           )
{
  p <- x/(x+a_inv)  ## in this manuscript p = 1-p
  if (sn1>0) tem <- p^y1*(1-p)^u1 else tem <- 1
  tem <- tem*dgamma(x, shape=sh, scale=sc)
  tem <- tem*pnbinom(y2-1, prob=1-p, size=u2)
  return(tem)
}



max.fun <-
  function(x=1,      #value of the random effect
           y1=0,                #sum(ypre); i.e. ypre+
           maxY2=3,             #max(Ynew)
           u1=1.5,              #Ex(Ypre+)
           u2=c(1.5,1.5,1.5),   #Ex(Ynew), vector
           a_inv=0.5,           #1/a
           sh=0.5, sc=2,sn1         #shame and scale of the Gamma dist'n
           )
{   
  P <- x/(x+a_inv) 
  ## pr(Y1+ = y1+ |g) 
  if (sn1 > 0) tem <- P^y1*(1-P)^u1 else tem <- 1
  tem <- tem*dgamma(x, shape=sh, scale=sc)
  for (i in 1:length(u2))
    { 
      tem1 <- pnbinom(maxY2-1, prob=1-P, size=u2[i])
      tem <- tem*tem1
    }
  return(tem) ## change from tem*v Jun 13
}



cum1.ln <-
function( x=2, #value of the random effect
         a_inv=0.5,             # 1/a
         sh=0.5, sc=2,          # lnorm paramters, sc=s=sqrt(log(Th+1)), sh = u = -s^2/2 
         y1=2, y2=2, u1=3, u2=3,sn1 #see cum.fun
         )
{   p=x/(x+a_inv)  
    if (sn1 > 0) tem <- p^y1*(1-p)^u1 else tem <- 1
    tem <- tem*dlnorm(x, meanlog=sh, sdlog=sc)
    tem <- tem*pnbinom(y2-1, prob=1-p, size=u2)
    return(tem)
}

max.ln <-
function(x=1, y1=0, maxY2=3, u1=1.5, u2=c(1.5,1.5,1.5), a_inv=0.5, 
       #see max.fun  
       sh=0.5, sc=2,sn1 #sh=mean, sc=sd of the lognormal dist'n
    )
{   
    P <- x/(x+a_inv) 
    #pr(Y1+= y1+ |g) 
    if (sn1 > 0) tem <- P^y1*(1-P)^u1 else tem <- 1
    tem <- tem*dlnorm(x, meanlog=sh, sdlog=sc)
    for (i in 1:length(u2)){
      
      tem1 <- pnbinom(maxY2-1, prob=1-P, size=u2[i])
      tem <- tem*tem1
    }
    return(tem)
}

## ===========================================

CP.se <-
  function(tpar,
           ## return
           ## 1) the point estimate of the conditional probability of observing
           ## the response counts as large as the observed ones given the previous counts
           ## 2) its asymptotic standard error
           
           ## estimates of log(alpha, sigma.G, betas); 
           ## e.g., tpar=olmeNB$est[,1] where olmeNB = output of the mle.*.fun function 
           Y1,          ## Y1 = y_i,old+ the sum of the response counts in pres
           Y2,          # Y2 = q(y_new)
           sn1, sn2,  # sn1 and sn2: number of scans in the pre and new sets.
           XM=NULL,       # XM : matrix of covariates
           dist="G",      # distribution of the random effects (G = gamma, N = log-normal)
           V,    # variance-covariance matrix of the parameter estimates; olmeNB$V
           pty="sum")     # q() = "sum" or "max" 
{ if (Y2==0) return(c(1,0))
  
  ## the point estimate Phat
  p <- jCP(tpar=tpar, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, XM=XM, dist=dist, type=pty, LG=FALSE)
  
  ## SE of logit(Phat)
  jac <- jacobian(func=jCP, x=tpar, Y1=Y1, Y2=Y2, sn1=sn1, sn2=sn2, XM=XM, dist=dist, LG=TRUE, type=pty)
  s2 <- jac%*%V%*%t(jac)
  s <- sqrt(s2) 
  return(c(p,s)) ## c(Phat, SE(logit(Phat)))
}


jCP <-
  function(tpar,
           ## estimates of log(alpha, sigma.G) betas; 
           ## e.g., tpar=olmeNB$est[,1] where olmeNB = output of the mle.*.fun function 
           Y1,
           ## Y1 = y_i,old+ the sum of the response counts in pres
           Y2,
           ## Y2 = q(y_new)
           sn1, sn2,
           ## sn1 and sn2: number of scans in the pre and new sets.
           XM=NULL,
           ## XM : matrix of covariates
           dist="G", 
           ## same as CP.se for description
           LG=FALSE,     #indicatior for logit transformation
           oth=NULL, # see Psum1 and Pmax1
           type="sum") 
{
  ## Return a point estimate
  a=exp(tpar[1])
  th=exp(tpar[2])

  sn=sn1+sn2 ## ni 
  
  u0=exp(tpar[3])
  u=rep(u0, sn)
  if (length(tpar)>3) u=u*exp(XM%*%tpar[-(1:3)])
  
  u1=sum(u[1:sn1])
  if (type=="sum") 
    { u2=sum(u[-(1:sn1)])
      temp=Psum1(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, Th=th, dist=dist, othr=oth,sn1=sn1)
    }else { ## max
      u2=u[-(1:sn1)]
      temp=Pmax1(Y1=Y1, Y2=Y2, u1=u1, u2=u2, a=a, Th=th, dist=dist, othr=oth,sn1=sn1)
    }
  if (LG) temp=lgt(temp)
  return(temp)
}


tp.fun <- function(i1, ## idat$ID
                   i2, ## idat$Vcode
                   y ## idat$CEL
                   )
{
  ## transpose the counts into a matrix so that each patient
  ## is a row and each visit is a column.
  ## It also puts an NA at any  visit with no data.
  x1 = table(i1,i2)
  x2 = match(x1,1) # 0->NA
  x = matrix(x2, nrow(x1), ncol(x1))
  dimnames(x) = dimnames(x1)
  i1 = as.numeric(factor(i1))
  i2 = as.numeric(factor(i2))
  for (i in 1:length(i1))
    {
      x[i1[i],i2[i]]=y[i]
      ## x1[i1[i],i2[i]]=prism.na$NASS.D[i]
    }
  ## return i1 by i2 contingency table
  return(x)
}


pp.ci <-
  function(x=c(0.01, 0.1), level=0.95, lg=TRUE)
  ## lg: logit transformation indicator
  ## x=c(phat, s.e.(phat)) if lg=F
  ## x=c(phat, s.e(logit(phat)) if lg=T
  ## level: confidence level
  ## x[1]: an estimate of conditional probability, p 
  ## x[2]: se (logit(phat)) if lg = TRUE
{
  if (is.na(x[1])){ return(c(NA, NA))
  }else if (x[1]==1) return(c(1,1))
                                        #lg=x[1]<0.2
  ll=0.5+level/2
  del=qnorm(c(1-ll, ll))
  tem=del*x[2]
  
  if (lg){
    return(ilgt(lgt(x[1])+tem))
  }else{
    tem1= x[1]+tem ##
    tem1[1] = max(tem1[1], 0)
    tem1[2] = min(tem1[2], 1)
    return(tem1)
  }
}





ppv.ci <- function(x = rbind(c(1, 0), c(0.1, 0.11)), level=0.95, lg=T)
                                        #see pp.ci
                                        # x is a n by 2 matrix
{   
  pp=cbind(x[,1], x[,1])
  ss=cbind(x[,2], x[,2])

  ll=0.5+level/2
  del=qnorm(c(1-ll, ll))
  tem=ss%*%diag(del)
  if (lg) 
    { r1 = pp*exp(tem)
      res=r1/(1-pp+r1)
    }else{ tem1= pp+tem 
           tem1[tem1[,1]<0,1] = 0
           tem1[tem1[,2]>1, 2] = 1
           res=tem1
         } #res[is.na(x[,1]), 1]=res[is.na(x[,1]), 2]= NA
                                        #res[x[,1]==1,2] = res[x[,1]==1,2] =1
  return(res)
}

pmarg.gauss.fun <-
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
    
    val=sum(pnbinom(y, size=uVa, prob=1-p)*GAU$weights)
    return(val)
  }

