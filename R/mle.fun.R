

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
           o.tol=1.e-3  # tolerance: for optim
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
  vcm <- solve(tt$hessian)   # the covariance matrix  = inverse Hessian
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


formulaToDat <-
function(formula, data, ID,labelnp=NULL) 
{

  datamatrix <- model.matrix(object=formula,data=data)
  covNames <- colnames(datamatrix)
  if ("(Intercept)" %in% colnames(datamatrix) )
    {
      datamatrix <- datamatrix[,-1,drop=FALSE]
    }else{
      stop("A model without an intercept term is not accepted!!")
    }
  y <- model.response(model.frame(formula=formula,data=data,na.action=na.omit))

  ## Update ID so that length(ID) = nrow(datamatrix) 
  formulaID <- update(formula,    ~ . + ID)
  dataID <- data; dataID$ID <- ID ## missing values in ID is not accepted! ID must be the same length as the n row of dataframe
  datamatrixID <- model.matrix(object=formulaID,data=dataID)
  upID <- datamatrixID[,colnames(datamatrixID)=="ID"] ## Updated ID does not have ID of missing observations
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
    return(list(dat=dat,labelnp=labelnp))
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


step.pl <-
function(gi, new=T, cl=1)
{ x=sort(unique(round(gi,4)))
  
  y=cumsum(table(gi))/length(gi)
  n=length(y)
  x=rep(x,rep(2,n))
  y=rep(y, rep(2,n))
  y=c(0,y)
  x=c(x,max(gi)+1)
  if (new)  plot(x,y, type="l", col=cl)
  else { lines(x,y, col=cl) }
  return()}


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
