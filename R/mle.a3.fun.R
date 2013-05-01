mle.a3.fun <-  function(formula,     ## an object of class "formula"
                        ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
                        data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
                        ## containing the variables in the model.
                        ID,          ## a vector of length n*ni containing patient IDs of observations in data
                        p.ini=NULL,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
                        IPRT=TRUE,      ## printing control
                        deps=1.e-4,   ## stop iteration when max(p.new - p.old) < deps
                        maxit = 100
                        )
{
  dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat)  
  #DT=list(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))

  if (is.null(p.ini))
  {
    p.ini=rep(0, 3+DT$cn) 
    p.ini[3]=mean(DT$y)   
  }
  p.new=p.ini
  W=NULL

  #initialize dth, ahat
  dth=1    #max(p.new-p.old)
  ahat=0
  counter <- 0 
  while(dth > deps & maxit > counter)
    {
    p.ini <- p.new           #reset the initial values at each iteration
    b.ini <- p.ini[-(1:2)]
    ## 
    ## Main Step 1 : estimate b using Generalized Least Square
    ## 
    if (DT$cn>0) #with covariates
      { temB <- optim(b.ini,## initial values of coefficients b0,b1,...
                      estb, ## weighted residuals function
                      hessian=T,
                      dat=DT, ## output of getDT
                      Wh=W, ## list of N elements, containing ni by ni weight matrix for each patient
                      PRT=F)
        bhat=temB$par
      }else{         #without covariates
        tem <- optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
        bhat <- tem$min
    }
    ## 
    ## Main Step 2: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)
    ## 
    
    ## Substep 1: the estimation of mu_i+
    ##         hat.mu_i+ = exp( X_i^T bhat)
    # compute exp(b0+b1*x1+ ...)
    uh=rep(exp(bhat[1]), DT$totN)
    if (DT$cn>0) {
      tem <- exp(DT$x%*%bhat[-1])
      uh <- tem*uh
    }
    #uhat=u.i+
    uhat <- as.vector(tapply(uh, DT$id, sum))

    ## Substep 2: the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
    #weights for gi
    if (ahat == 0) #first iteration, set wi=1
    {
      wi <- rep(1, DT$np)
    }else{
      gsig <- that+(1+ahat*(that+1))/uhat
      wi <- sqrt(that/gsig)   
    }
    
    gh0 <- DT$ys/uhat ## An estimate of y_i+/mu_i+ for each i
 
    gh=wi*gh0+1-wi ## An estimate of w_i*(y_i+/mu_i+)+(1-w_i) for each i
    
    if (IPRT) { ##print(summary(wi)) ; print(summary(gh0));
      cat("\n iteration:",counter)
      cat("\n The distribution of hat.g \n")
      print(summary(gh))
    }
    
    #normalized gh so that the average is 1
    gh=gh/mean(gh)
    gh=as.vector(gh)

    # frequency table for gh
    tem <- sort(round(gh,6))
    gh1 <- unique(tem)
    ghw <- table(tem)/DT$np
    ## the number of unique values of gh by 1 vector, containing the proportion of observations with that value
    gtb=cbind(ghw, gh1) 
    rownames(gtb)=NULL
    ## 
    ## Main Step 3: estimate a using profile likelihood
    ## 
    tt <- optimize(lk.a.fun, lower=-5, upper=5, dat=DT, Iprt=F, gs=gtb,
                   uhat=uh ## a vector of length n*sn containing hat.mu_ij 
                   )
    ## lk.a.fun is a function of log(a)
    ## exponentiate back 
    ahat=exp(tt$min)
    ## 
    ## Main Step 4: moment estimate of var(G)
    ##
    that=estSg3(DT, uhat, uh) ## uh is a vector of length length(Y)

    #weights for GLS
    W=tapply(uh, DT$id, getWhs, th=that, a=ahat)
    
    p.new=c(tt$min, log(that), bhat)
    if (IPRT) {cat("\nlog(a), log(th), b0, b1, ...\n");print(p.new)}
    dth=max(abs(p.new-p.ini))
    counter <- counter +1
  } 
  vcm=NULL 
  if (DT$cn>0) vcm=solve(temB$hessian)
  if (counter == maxit) warning("The maximum number of iterations occur!")
  p.new <- matrix(p.new,ncol=1)
  rownames(p.new) = c("log_a", "log_var(G)", "(Intercept)", DT$xnames)
  if (is.matrix(vcm)) colnames(vcm) <- rownames(vcm) <- c("(Intercept)", DT$xnames)
  re <- list(opt = tt, V = vcm, est = p.new, gtb = gtb,
              gi = as.vector(gh), mod = "NoN", 
              ##idat = data.frame(dat),
              cor="ind",formula=formula)
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}


estb <-
  function(b, # (b0, b1, ...)
           dat, #data see mle.fun
           Wh=NULL,
           ## a list of N elements.
           ## each element is ni by ni matrix of weights
           ## If Wh=NULL then the weight matrix are all regarded as identity
           PRT=T)
{
  ## GLS compute the weighted sum of the residuals
  th2=exp(b[1]) #th2 =exp(b0)
  ## If the number of covariate is nonzero:
  if (dat$cn>0) 
  {
    tem=exp(dat$x%*%b[-1])
    th2=tem*th2 #th2 = exp(b0+b1*x1 + ...) n*sn by 1 vector
  }
  tem=(dat$y-th2)
  ## If weights are not provided then  W_i = diag(dt$vn) identity
  ## so the weighted sum of residuals is
  ## sum_i (y_i-X_i^T beta)^T (y_i-X_i^T beta) = sum(tem^2)
  if (is.null(Wh))
  { 
    val=sum(tem^2)
  }
  else
  { val=0
    for (i in 1:dat$np)
    {
      ## ii is a vector of length ni containing the locations of
      ## repeated measures of patient i
      ii=dat$ind[i]:(dat$ind[i+1]-1)
      val=val+t(tem[ii])%*%Wh[[i]]%*%tem[ii]
    }
  }
  if (PRT) print(val)
  return(val)
}

estSg3 <-
function(
         dat,  ## data from getDT
         ui,   ## u.i+
         uij   ## u.ij, vector
         )
{ 
  ssy=sum(dat$ys^2)-sum(dat$y^2) ## SS
  ssu=sum(ui^2)-sum(uij^2)
  that=ssy/ssu
  return(max(that-1,0)) 
}


lk.a.fun <-
  function(para=-0.5, #log(a)
           dat, Iprt=FALSE, 
           gs, #gi in freqency table form
           uhat    #hat u.ij, same length as Y.ij
           ) 
{
  ## psudo profile likelihood as a function of alpha
  ## Li(a; hatb, hatg, y)
  ## = 1/N sum_l=1^N prod_j=1^ni Pr(Y_ij=y_ij|G_l=g_l,hatb,a)
  ## = sum_l=1^N* prod_j=1^ni Pr(Y_ij=y_ij|G_l=g*_l,hatb,a)*prop(g*_l)
  ## = sum_l=1^N* prod_j=1^ni choose(y_ij+r_ij-1,y_ij) (1-p(g*_l,a))^r_ij p(g*_l,a)^{y_ij}*prop(g*_l)
  ## = [prod_j=1^ni choose(y_ij+r_ij-1,y_ij)] sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)
  ## 
  ## log Li
  ## = [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] + log [sum_l=1^N* (1-p(g*_l,a))^r_i+ p(g*_l,a)^{y_i+}*prop(g*_l)]
  ## 
  ## g*_l represents the unique value of the approximations
  ## prop(g*_l) represents the proportion of approximated g_l's that have the value g*_l
  ## N* represents the number of unique g*
  ## p = p(g*_l,a) = 1/(g*_l*a + 1)
  
  if (Iprt) cat(para)
  
  a=exp(para)
  
  ainv=1/a  
  
  th2=uhat/a ## r_ij = exp(X_ij^T hat.beta)/alpha

  ## -log [sum_j=1^ni log choose(y_ij+r_ij-1,y_ij)] = sum( - lgamma(dat$y+th2) + lgamma(th2) + lgamma(dat$y))
  ##  but the last term is not necessary because it does not depend on a
  nllk=sum( - lgamma(dat$y+th2) + lgamma(th2))
  us=tapply(th2, dat$id, sum)
  
  p=gs[,2]/(gs[,2]+ainv)
  ## please note that p = 1 - p in this manuscript
  ## gs[,2] is a vector of length # unique g_i containing the values of unique g_i
  
  q=1-p
  for (i in 1:dat$np)
    {  tem=sum(
         p^dat$ys[i]*q^us[i]*gs[,1]
         ## gs[,1] is a vector of length # unique g_i
         ## containing the proportion of g_i
         )
      nllk=nllk-log(tem)
     }
  
  if (Iprt) cat(" nllk=", nllk, "\n")
  return(nllk)
}


getWhs <-
  function(
           u=c(1.5, 1.35, 1.2, 1.1, 1, 0.9), #u.ij, vector
           th=exp(1.3), #var(G)
           a=exp(-0.5)  #a
           )
{
  ## Independent model has:
  ## Var(Y_ij)=mu_ij^2*var(G) + mu_ij*(1+(var(G)+1)*a)
  w1=u%*%t(u)*th
  if (length(u)==1)w2=(1+(th+1)*a)*u ## Modified by Yumi 2013 Mar8
    else w2=(1+(th+1)*a)*diag(u)
  W=solve(w1+w2)
  return(W)
}



mle.a1.fun <-
function(dat, p.ini=NULL, IPRT=TRUE, deps=1e-4)
 { 
  DT=getDT(dat)  
  #DT=(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))

  if (is.null(p.ini))
  { p.ini=rep(0, 3+DT$cn) 
    p.ini[3]=mean(DT$y)
  }
  p.new=p.ini
  W=NULL

  dth=1 #max(p.new-p.ini)
  while(dth > deps)
  { p.ini=p.new
    b.ini=p.ini[-(1:2)]
    
    if (DT$cn>0) 
    { temB=optim(b.ini, estb, hessian=T, dat=DT, Wh=W, PRT=F)
      bhat=temB$par }
    else
    { tem=optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
      bhat=tem$min
    }

    uh=rep(exp(bhat[1]), DT$totN)
    if (DT$cn>0) { tem=exp(DT$x%*%bhat[-1])
    uh=tem*uh }

    uhat=as.vector(tapply(uh, DT$id, sum))
    gh0=DT$ys/uhat

    #for cases of gh0=0, uniformly between 0 to 1/2 of the smallest non zero value
    gh_1=min(gh0[gh0>0])/2

    ll=sum(gh0==0)
 
    gh=gh0
    ll1=ll%%10
    ll2=ll%/%10
    ll3=rep(1:0, c(ll1, 10-ll1))+ll2
    gh[gh0==0] = rep(seq(0, gh_1, length=10), ll3)
   
    gh=gh/mean(gh)
  
    gh=as.vector(gh)

    #frequency table
    tem=sort(round(gh,6))
    gh1=unique(tem)
    ghw=table(tem)/DT$np
    gtb=cbind(ghw, gh1)
    rownames(gtb)=NULL

    tt=optimize(lk.a.fun, lower=-5, upper=5, dat=DT, Iprt=F, gs=gtb, uhat=uh)
    ahat=exp(tt$min)

     that=estSg3(DT, uhat, uh)
     W=tapply(uh, DT$id, getWhs, th=that, a=ahat)

    p.new=c(tt$min, log(that), bhat)
    if (IPRT) print(p.new)
    dth=max(abs(p.new-p.ini))
  }  

  if (DT$cn>0) vcm=solve(temB$hessian)
  names(p.new)=c("log_a", "log_th", "log_u0", DT$xnames)
  return(list(opt=tt, vcm=vcm, est=p.new, gi=as.vector(gh), gtb=gtb, mod="NoN"))
}


## mle.a3.fun <-
##   function(formula,     ## an object of class "formula"
##            ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##            data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
##            ## containing the variables in the model.
##            ID,          ## a vector of length n*ni containing patient IDs of observations in data
##            p.ini=NULL,  ## initial values for the parameters (log(a), log(th), b0, b1, ...)
##            IPRT=TRUE,      ## printing control
##            deps=1.e-4   ## stop iteration when max(p.new - p.old) < deps
##            )
## {
##   dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
##   DT=getDT(dat)  
##   #DT=list(id, y, x, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND, xnames=xnames))

##   if (is.null(p.ini))
##   {
##     p.ini=rep(0, 3+DT$cn) 
##     p.ini[3]=mean(DT$y)   
##   }
##   p.new=p.ini
##   W=NULL

##   #initialize dth, ahat
##   dth=1    #max(p.new-p.old)
##   ahat=0
##   counter <- 0 
##   while(dth > deps)
##     {
##     p.ini=p.new           #reset the initial values at each iteration
##     b.ini=p.ini[-(1:2)]
##     ## 
##     ## Main Step 1 : estimate b using Generalized Least Square
##     ## 
##     if (DT$cn>0) #with covariates
##       { temB=optim(b.ini,## initial values of coefficients b0,b1,...
##           estb, ## weighted residuals function
##           hessian=T,
##           dat=DT, ## output of getDT
##           Wh=W, ## list of N elements, containing ni by ni weight matrix for each patient
##           PRT=F)
##         bhat=temB$par
##       }else{         #without covariates
##         tem=optimize(estb, lower=-5, upper=5, dat=DT, Wh=W, PRT=F)
##         bhat=tem$min
##     }
##     ## 
##     ## Main Step 2: estimate g_i i=1,...,N, g_i = w_i*(y_i+/mu_i+)+(1-w_i)
##     ## 
    
##     ## Substep 1: the estimation of mu_i+
##     ##         hat.mu_i+ = exp( X_i^T bhat)
##     # compute exp(b0+b1*x1+ ...)
##     uh=rep(exp(bhat[1]), DT$totN)
##     if (DT$cn>0) { tem=exp(DT$x%*%bhat[-1])
##     uh=tem*uh }
##     #uhat=u.i+
##     uhat=as.vector(tapply(uh, DT$id, sum))

##     ## Substep 2: the estimation of w_i = sqrt(var(G_i)/var(Y_i+/mu_i+))
##     #weights for gi
##     if (ahat == 0) #first iteration, set wi=1
##     {
##       wi=rep(1, DT$np)
##     }else{
##       gsig=that+(1+ahat*(that+1))/uhat
##       wi=sqrt(that/gsig)   
##     }
    
##     gh0=DT$ys/uhat ## An estimate of y_i+/mu_i+ for each i
 
##     gh=wi*gh0+1-wi ## An estimate of w_i*(y_i+/mu_i+)+(1-w_i) for each i
    
##     if (IPRT==T) { ##print(summary(wi)) ; print(summary(gh0));
##       cat("\n iteration:",counter)
##       cat("\n The distribution of hat.g \n")
##       print(summary(gh))
##     }
    
##     #normalized gh so that the average is 1
##     gh=gh/mean(gh)
##     gh=as.vector(gh)

##     # frequency table for gh
##     tem=sort(round(gh,6))
##     gh1=unique(tem)
##     ghw=table(tem)/DT$np
##     ## the number of unique values of gh by 1 vector, containing the proportion of observations with that value
##     gtb=cbind(ghw, gh1) 
##     rownames(gtb)=NULL
##     ## 
##     ## Main Step 3: estimate a using profile likelihood
##     ## 
##     tt=optimize(lk.a.fun, lower=-5, upper=5, dat=DT, Iprt=F, gs=gtb,
##       uhat=uh ## a vector of length n*sn containing hat.mu_ij 
##       )
##     ## lk.a.fun is a function of log(a)
##     ## exponentiate back 
##     ahat=exp(tt$min)
##     ## 
##     ## Main Step 4: moment estimate of var(G)
##     ##
##     that=estSg3(DT, uhat, uh) ## uh is a vector of length length(Y)

##     #weights for GLS
##     W=tapply(uh, DT$id, getWhs, th=that, a=ahat)
    
##     p.new=c(tt$min, log(that), bhat)
##     if (IPRT) {cat("\nlog(a), log(th), b0, b1, ...\n");print(p.new)}
##     dth=max(abs(p.new-p.ini))
##     counter <- counter +1
##   } 
##    vcm=NULL 
##    if (DT$cn>0) vcm=solve(temB$hessian)
  
##   names(p.new) = c("log_a", "log_var(G)", "beta0", DT$xnames)

##   return(list(opt = tt, vcm = vcm, est = p.new, gtb = gtb,
##               gi = as.vector(gh), mod = "NoN", 
##               idat = data.frame(dat),cor="ind"))
## }

