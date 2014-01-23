

lmeNB <- function(formula,data,ID,p.ini=NULL,IPRT=FALSE,AR=FALSE,RE=c("G","N","semipara"),deps=1e-03,Vcode,
                  i.tol=1e-75,
                  o.tol=1.e-3,
                  maxit=100)
  {
    ## A Wrapper function to call
    ## mle.fun
    ## mle.a3.fun
    ## mle.ar1.fun
    ## mle.ar1.non3
    if (length(RE)>1) RE <- RE[1]
    if (!AR){
      
      if (RE=="G" || RE == "N"){
        return( mle.fun(formula=formula,data=data,ID=ID,p.ini=p.ini,IPRT=IPRT,model=RE))
      }else if (RE=="semipara"){
        return(mle.a3.fun(formula=formula,data=data,ID=ID, p.ini = p.ini, IPRT = IPRT, deps = deps,maxit=maxit))
      }else{
        stop("RE must be G, N or semipara!!")
      }
      
    }else{
      
      if (RE=="G" || RE == "N"){
        return(mle.ar1.fun(formula=formula, data=data, ID=ID, Vcode=Vcode,
                           p.ini=p.ini, IPRT = IPRT, model = RE, 
                           i.tol = i.tol, o.tol = o.tol
                            ))
      }else if (RE=="semipara"){
        return(mle.ar1.non3(formula=formula, data=data, ID=ID, Vcode=Vcode,
                            p.ini = p.ini, IPRT = IPRT, deps = deps, maxit=maxit))
      }else{
        stop("RE must be G, N or semipara!!")
      }
      
    }
  }


mle.fun <-
  function(formula,     ## an object of class "formula"
           ## (or one that can be coerced to that class): a symbolic description of the model to be fitted.
           data,        ## a data frame, list or environment (or object coercible by as.data.frame to a data frame)
           ## containing the variables in the model.
           ID,          ## a vector of length n*ni containing patient IDs of observations in data
           p.ini=NULL,    # initial value c(log(a), log(th), b0, b1, ...)
           IPRT=FALSE,   # printing control: T = print iterations
           model="G", #RE options: G = gamma; N = lognormal
           i.tol=1e-75,
           o.tol=1.e-3,  # tolerance: for optim
           COV=TRUE     ## Covariance matrix == TRUE
           )
{ ##########################
  ## Tips for computation ##
  ##########################
  ## Note that mle.fun fails to compute the negative log-likelihood values
  ## when the total sum of response counts of a patient is **very large**
  ## This is because the log-likelihood is:
  ## 
  ##         integrate dnbinom(sum_j=1^n_i y_ij; sum_j=1^n_i r_ij,p_i)*f(g) dg
  ## \propto integrate p_i^{sum_j=1^n_i r_ij} (1-p_i)^{sum_j=1^n_i y_ij} *f(g) dg
  ##
  ## since 0 <1-p_i < 1 if sum_j=1^n_i y_ij is **very** large, then (1-p_i)^{sum_j=1^n_i y_ij} = 0.
  ## and the integrated value become zero. e.g., 0.5^1000 = 0 in R
  ## Since log of this quantity is taken,
  ## the evaluated value become log(0) = - Inf
  ##
  cmd <- match.call() 
  dat <- formulaToDat(formula=formula,data=data,ID=ID)  # dat = (ID, Y, x1, x2, ...) numeric matrix
  DT <- getDT(dat)

  if (is.null(p.ini))
  {
    p.ini=rep(0, 3+DT$cn) 
    p.ini[3]=mean(DT$y)   
  }
  
  if (IPRT) cat("Print log_a log_th log_u0", DT$xnames, " and negative of the log-likelihood at each iteration")
  
  if (model=="G") ## Gi ~ Gamma(scale=theta,shape=1/theta)
   tt <- optim(p.ini, lk.fun,control=list(reltol=o.tol), hessian=TRUE, dat=DT, Iprt=IPRT,tol=i.tol)
  else if (model=="N")## Gi ~ logN(mu=1,var(G_i)=theta)
  tt <- optim(p.ini, lk.ln.fun,control=list(reltol=o.tol), hessian=TRUE, dat=DT, Iprt=IPRT,tol=i.tol)
    
  nlk <- tt$value+sum(lgamma(DT$y+1))  #full -loglikelihood
  ## if (IPRT) print(tt$hessian)
  if (COV){
    vcm <- solve(tt$hessian)   # the covariance matrix  = inverse Hessian
  }else vcm <- matrix(NA,length(tt$p),length(tt$p))
  p.est <- cbind(tt$p, sqrt(diag(vcm)))
  row.names(p.est) <- rownames(vcm)<-colnames(vcm)<- c("log_a", "log_th", "(Intercept)", DT$xnames)
  colnames(p.est) <- c("estimate","SE")
  
  re <- list(call=cmd, p.ini=p.ini,opt=tt,formula=formula,
             nlk=nlk, V=vcm, est=p.est, mod=model, ##idat=data.frame(dat),
             cor="ind")
  class(re) <- "LinearMixedEffectNBFreq"
  return(re)
}

print.LinearMixedEffectNBFreq <- function(x,...)
  {

    cat("\n -----Negative binomial mixed effect regression----- ")
    if (x$mod == "G") RE <- "gamma"
    else if (x$mod == "N") RE <- "log-normal"
    else if (x$mod == "NoN") RE <- "semi-parametric"
    if (x$cor=="ind") depn <- "independent"
    else if (x$cor=="ar1") depn <- "AR(1)"
    cat("\n ---Random effect is from ",RE," distribution---")
    cat("\n The ",depn," correlation structure for the conditional dist'n of response given random effects.")

    cat("\n Formula: ");print(x$formula)
    cat("\n Estimates: \n")
    crit <- qnorm(0.975)
    if (x$mod != "NoN")
      {
        printForm <- data.frame(cbind(sprintf("%1.3f",x$est[,1]),sprintf("%1.3f",x$est[,2]),
                                     sprintf("%1.3f",x$est[,1]-crit*x$est[,2]),
                                     sprintf("%1.3f",x$est[,1]+crit*x$est[,2])))
        rownames(printForm) <-rownames(x$est)
        colnames(printForm) <-c("Value","Std.Error","lower CI","upper CI")
      }else{
        printForm <- x$est
        colnames(printForm) <- "Value"
      }
    print(printForm) 
    cat("\n Estimated covariance matrix: \n")
    print(x$V)
    cat("------------------------")
    if (x$mod != "NoN") cat("\n Log-likelihood",-x$nlk)
  }


formulaToDat <- function(formula, data, ID,labelnp=NULL) 
{
  ## datamatrix exclude the covariate values corresponding to the NA CELs or NA selected covariates
  datamatrix <- model.matrix(object=formula,data=data)
  covNames <- colnames(datamatrix)
  if ("(Intercept)" %in% colnames(datamatrix) )
    {
      datamatrix <- datamatrix[,-1,drop=FALSE]
    }else{
      stop("A model without an intercept term is not accepted!!")
    }
  ## y exclude the covariate values corresponding to the NA CELs or NA selected covariates
  y <- model.response(model.frame(formula=formula,data=data,na.action=na.omit))

  ## Update ID so that length(ID) = nrow(datamatrix)
  formulaNonMiss <- update(formula,    ~ . + NonMiss)
  dataNonMiss <- data; dataNonMiss$NonMiss <- 1:nrow(data)
  datamatNonMiss <- model.matrix(object=formulaNonMiss,data=dataNonMiss)
  upID <- origID <- ID[ 1:nrow(data) %in% datamatNonMiss[,colnames(datamatNonMiss)=="NonMiss"] ]
  
  ## If ID is a character vector of length sum ni,
  ## it is modified to an integer vector, indicating the first appearing patient
  ## as 1, the second one as 2, and so on..
   
  uniID <- unique(upID)
  ID <- rep(NA,length(upID))
  for (i in 1 : length(uniID))
    {
      ID[upID == uniID[i]] <- i

      s <- sum(upID == uniID[i])
      for (j in 1 : s )
        upID[upID == uniID[i]] <- paste(uniID[i],1:s,sep="--")
    }
  dat <- cbind(ID=ID,CEL=y,datamatrix)
  rownames(dat) <- upID

  if (! is.null(labelnp)){
    formulalnp <- update(formula, ~.+labelnp)
    datalnp <- data
    datalnp$labelnp <- labelnp
    mm <- model.matrix(object=formulalnp,data=datalnp)
    labelnp <- mm[,colnames(mm)=="labelnp"]
    return(list(dat=dat,labelnp=labelnp,origID=origID))
  }
  return(dat) #dat=(ID, Y, x1, x2, ...)
}

getDT <- function(dat)             #dat=(ID, Y, x1, x2, ...)
{
  dat <- dat[order(dat[,1]),]
  ID <- dat[,1]
  ydat <- dat[,2]
  ncov <- ncol(dat)-2
   
  if (ncov>0){
    xdat=as.matrix(dat[,3:ncol(dat)])
    xnames=colnames(dat)[-(1:2)] 
    ##print(xnames)
  }else{
    xdat=NULL
    xnames=NULL
  } 

  Ysum <- tapply(ydat, ID, sum) #total lesion count for each patient

  Ni <- table(ID) # number of scans for each patient
  totN <- sum(Ni) # total number scans of cohort
  IND <- c(0,cumsum(Ni))+1 #location of the 1st row for each patient
  N <- length(Ni) #number of patients
  
  return(list(id=ID, y=ydat, x=xdat, cn=ncov, ys=Ysum, np=N, totN=totN, ni=Ni,  ind=IND, xnames=xnames))
}


ilgt <-
function(x)
{tem=exp(x)
 res=tem/(1+tem)
 return(res)
 }


## Likelihood function for gamma random effect model
lk.fun <-
  function(para, ## parameters (log(a), log(th), b0, b1, ...) 
           dat,  ## dat=list(id=ID, y=ydat, x=xdat, cn=ncov, 
           ##                ys=Ysum, np=N, totN=totN, ni=Ni, ind=IND)
           ## output of getDT()
           Iprt=FALSE, #printing control 
           tol=1.e-75, #control parameter for integration 
           sig=FALSE # if T, compute the full likelihood 
                                        #(i.e., with the constant)
           )
{ if (Iprt) cat(para)
  a <- exp(para[1])
  th1 <- exp(para[2]) #scale
  ainv <- 1/a  
  shp <- 1/th1

  u0 <- exp(para[3])
  th2 <- rep(u0/a, dat$totN)
  if (dat$cn>0) {
    ## if there are covariates
    b <- para[4:(dat$cn+3)]
    tem <- exp(dat$x%*%b)
    th2 <- tem*th2 ## th2 = r_{ij}=exp(X^T beta)/a: Y_{ij}|G_i=gi ~NB(r_{ij},p_i)
  }
  ## -loglikelihood
  ## constant terms of - loglikelihood
  nllk <- sum( - lgamma(dat$y+th2) + lgamma(th2))
  us <- tapply(th2, dat$id, sum)
  lk <- rep(0, dat$np)
  for (i in 1:dat$np)
    { tem=integrate(int.fun, lower=0, upper=Inf, a_inv=ainv, abs.tol=tol,
       ## dat$ys[i] = sum(y_{ij},j=1,...,ni)
            sh=shp, sc=th1, ysum=dat$ys[i],  usum=us[i])
     #print(c(dat$ys[i], us[i], tem$v))
     nllk = nllk - log(tem$v)
     #sig=T full likelihood
     if (sig) 
     {  ll=dat$ind[i]:(dat$ind[i+1]-1)
        lk[i]=log(tem$v)+sum(lgamma(dat$y[ll]+th2[ll]) - lgamma(th2[ll])-lgamma(dat$y[ll]+1))
     }
   }
   if (Iprt) cat(" nllk=", nllk, "\n")
   if (sig) return(lk) #return full likelihood (vector of length n)
   else return(nllk)   #return log likelihood without constant
}

int.fun <-
function(x=2, # value of G 
         a_inv=0.5,  # a_inv=1/a
         sh=0.5, sc=2, #shape and scale; sh=1/th; sc=th to have E(G_i)=scale*shape=1
         ysum=2, #sum(y.ij,j=1,..,ni); ysum is a number
         usum=3  #sum(u.ij/a,j=1,...,ni); usum is a number
         )
  ## note p=1-p
{   p=x/(x+a_inv)  
    tem=p^ysum*(1-p)^usum*dgamma(x, shape=sh, scale=sc)
    return(tem)
}


## step.pl <-
## function(gi, new=T, cl=1)
## { x=sort(unique(round(gi,4)))
  
##   y=cumsum(table(gi))/length(gi)
##   n=length(y)
##   x=rep(x,rep(2,n))
##   y=rep(y, rep(2,n))
##   y=c(0,y)
##   x=c(x,max(gi)+1)
##   if (new)  plot(x,y, type="l", col=cl)
##   else { lines(x,y, col=cl) }
##   return()}


## ========= log-normal random effects =========
lk.ln.fun <-
function(para, dat, Iprt=F, tol=1e-75)
{ if (Iprt) cat("\n",para)
  
  a=exp(para[1])
  th1=exp(para[2]) #scale
  
  s.ln=log(th1+1) #s^2
  u.ln=-s.ln/2    #-s^2/2, mean for dlnorm()
  s.ln=sqrt(s.ln) #s, sd for dlnorm()

  ainv=1/a  

  u0=exp(para[3])
  th2=rep(u0/a, dat$totN)
  if (dat$cn>0) {
   b=para[4:(dat$cn+3)]
   tem=exp(dat$x%*%b)
   th2=tem*th2
  }

   nllk=sum( - lgamma(dat$y+th2) + lgamma(th2))
   us=tapply(th2, dat$id, sum)
   
   for (i in 1:dat$np) 
   { tem=integrate(int1.ln, lower=0, upper=Inf, a_inv=ainv,
            sh=u.ln, sc=s.ln, ysum=dat$ys[i],  usum=us[i], abs.tol=tol)
     #print(c(dat$ys[i], us[i], tem$v))
     nllk = nllk - log(tem$v)
   }
    if (Iprt) cat(" nllk=", nllk)
   return(nllk)
}




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

    ## If p.new and p.ini contains Inf at the same spot, then
    ## dth = NaN
    whereInf <- (abs(p.new)==Inf)* (abs(p.new)==Inf)
    if ( sum(whereInf))
      {
        dth=max(abs(p.new[whereInf]-p.ini[whereInf]))
      }
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

