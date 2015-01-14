\name{fitParaAR1}
\alias{fitParaAR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs the maximum likelihood estimation for the negative binomial mixed-effect AR(1) model 
}

%
\description{
This function fits a negative binomial mixed-effect AR(1) model in the formulation described Zhao et al. (2013).
The conditional distribution of response counts given random effect is modelled by Negative Binomial as described in description of \code{\link{lmeNB}}. 
The conditional dependence among the response counts of a subject is modeled with AR(1) structure. 
The random effects are modelled with either gamma or log-normal distributions. 
See descriptions of \code{\link{lmeNB}}.
}

\usage{
fitParaAR1(formula, data, ID, Vcode, p.ini = NULL, IPRT = FALSE, 
           RE = "G", i.tol = 1e-75, o.tol = 0.001) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{


\item{formula}{
See \code{\link{lmeNB}}.
}
\item{data}{
See \code{\link{lmeNB}}.
}
\item{ID}{
See \code{\link{lmeNB}}.
}
\item{Vcode}{
See \code{\link{lmeNB}}.
}

  \item{p.ini}{
A vector of length 4 + # covariates, containing the initial values for the parameters 
 (\eqn{\log(\alpha)}{log(\alpha)},
\eqn{\log(\theta)}{log(\theta)},
\eqn{logit(\delta)}{logit(\delta)},
\eqn{\beta_0}{\beta[0]}, 
\eqn{\beta_1...}{\beta[1]...}).
\code{NULL} is accepted.
}

\item{IPRT}{	
See \code{\link{lmeNB}}.
}
  \item{RE}{
See \code{\link{fitParaIND}}.
}

\item{i.tol}{
See \code{\link{lmeNB}}.
}
  \item{o.tol}{
See \code{\link{lmeNB}}.
}
  

}

%to add
\details{
\code{fitParaAR1}  calls \code{optim} to minimize the negative log-likelihood of the negative binomial model with respect to the model parameters:
(\eqn{\log(\alpha)}{log(\alpha)},
\eqn{\log(\theta)}{log(\theta)},
\eqn{logit(\delta)}{logit(\delta)},
\eqn{\beta_0}{\beta[0]}, 
\eqn{\beta_1...}{\beta[1]...}).
The Nelder-Mead algorithm is employed.
The log-likelihood is obtained by marginalizing out the random effects. 
The numerical integration is carried out using adaptive quadrature.
When missing visits are present, an approximation of the likelihood is used (see Zhao et al. (2013) for details.)
All the computations are done in \code{R}.
}

\value{

\item{opt}{
See \code{\link{lmeNB}}.
}
\item{nlk}{
See \code{\link{lmeNB}}.
}
\item{V}{ 
See \code{\link{lmeNB}}.
}
\item{est}{ 
See \code{\link{lmeNB}}.
}
%%\item{RE}{}
%%\item{Vcode}{}
\item{AR}{\code{TRUE}}
%%\item{formula}{}
}

\references{
Detection of unusual increases in MRI lesion counts in individual multiple sclerosis patients. (2013) Zhao, Y., Li, D.K.B., Petkau, A.J., Riddehough, A., Traboulsee, A., Journal of the American Statistical Association.
}
\author{
 Zhao, Y. and Kondo, Y.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The main function to fit the Negative Binomial mixed-effect model:
\code{\link{lmeNB}},

The functions to fit the other models:
\code{\link{fitParaIND}},
%\code{\link{fitParaAR1}},
\code{\link{fitSemiIND}},
\code{\link{fitSemiAR1}},

The subroutines of \code{\link{index.batch}} to compute the conditional probability index:
\code{\link{jCP.ar1}},
\code{\link{CP1.ar1}},
\code{\link{MCCP.ar1}},
\code{\link{CP.ar1.se}},
\code{\link{CP.se}},
\code{\link{jCP}},

The functions to generate simulated datasets:
\code{\link{rNBME.R}}.
}

\examples{
\dontrun{

## ==================================================================================
## generate a data based on the negative binomial mixed-effect AR(1) model.
## Under this model, the response counts follows the negative binomial:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ Gamma(scale=th,shape=1/th)
## 
## The adjacent repeated measures of the same subject are correlated 
## with correlation structure:
## cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)  

loga <- -0.5
logtheta<- 1.3
logitd <- -0.2
b0 <- 0.5 ## no covariates; 
## 80 subjects each with 5 scans
n <- 80
sn <- 5

set.seed(1)
DT2 <-  rNBME.R(gdist = "G",
               n = n, ## 	the total number of subjectss	       
	       sn = sn,
	       th=exp(logtheta),
               u1 = rep(exp(b0),sn),
	       u2 = rep(exp(b0),sn),
	       a = exp(loga),
	       d = exp(logitd)/(1+exp(logitd))
	      )
Vcode <- rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
ID <- DT2$id
new <- Vcode > 0
dt2 <- data.frame(CEL=DT2$y)

## ================================================================================

## 1) Fit the negative binomial mixed-effect AR(1) model 
## where the random effects are from the gamma distribution
## This is the true model

re.gamma.ar1 <- fitParaAR1(formula=CEL~1,data=dt2,ID=ID,
		           Vcode=Vcode, 
		           p.ini=c(loga,logtheta,logitd,b0), 
		           ## log(a), log(theta), logit(d), b0
		           RE="G", 
		           IPRT=TRUE) 



## compute the estimates of the conditional probabilities 
## with sum of the new repeated measure as a summary statistics 
## Note C=TRUE with i.tol=1E-3 options compute the index faster
## i.se=TRUE requires more time
Psum <- index.batch(olmeNB=re.gamma.ar1,data=dt2,ID=ID,Vcode=Vcode,
	 	   labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE,C=TRUE,i.tol=1E-3)  
		 

## compute the estimates of the conditional probabilities 
## with max of the new repeated measure as a summary statistics 
Pmax <-index.batch(olmeNB=re.gamma.ar1,data=dt2,ID=ID,Vcode=Vcode, 
                  labelnp=new,qfun="max", IPRT=TRUE,i.se=FALSE,C=TRUE,i.tol=1E-3)

## Which patient's estimated probabilities based on the sum and max 
## statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt2$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,]


## 2) Fit the negative binomial mixed-effect AR(1) model 
## where random effects is from the log-normal distribution

re.logn.ar1 <- fitParaAR1(formula=CEL~1,data=dt2,ID=ID,
		          Vcode=Vcode, 
		          p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
   		          RE="N",IPRT=TRUE)

Psum <- index.batch(olmeNB=re.logn.ar1,data=dt2,ID=ID,Vcode=Vcode, 
                    labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE,C=TRUE,i.tol=1E-3) 
re.logn.ar1$Psum <- Psum


## 3) Fit the negative binomial independent model 
## where random effects are from the gamma distribution
re.gamma.ind <- fitParaIND(formula=CEL~1,data=dt2,ID=ID, 
                           RE="G", 
	         	   p.ini=c(loga,logtheta,b0), 
		           IPRT=TRUE)

Psum <- index.batch(olmeNB=re.gamma.ind,data=dt2,ID=ID, 
                    labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE)  



## 4) Fit the negative binomial independent model 
## where random effects are from the lognormal distribution
re.logn.ind <- fitParaIND(formula=CEL~1,data=dt2,ID=ID, 
                          RE="N", 			   	
		          p.ini=c(loga,logtheta,b0), 		
		          IPRT=TRUE)

Psum <- index.batch(olmeNB=re.logn.ind, data=dt2,ID=ID, 
                   labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE) 


## 5) Fit the semi-parametric negative binomial AR(1) model 

logvarG <- -0.5
re.semi.ar1 <- fitSemiAR1(formula=CEL~1,data=dt2,ID=ID, 
                          p.ini=c(loga, logvarG, logitd,b0),Vcode=Vcode)
Psum <- index.batch(olmeNB=re.semi.ar1,data=dt2,ID=ID, Vcode=Vcode,
	       	    labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE)  


## 6) Fit the semi-parametric negative binomial independent model 
re.semi.ind <- fitSemiIND(formula=CEL~1,data=dt2,ID=ID, p.ini=c(loga, logvarG, b0))
Psum <- index.batch(olmeNB=re.semi.ind,data=dt2,ID=ID,  
                   labelnp=new, qfun="sum", IPRT=TRUE,i.se=FALSE) 


## ======================================================================== ##
## == Which model performed the best in terms of the estimation of beta0 == ##
## ======================================================================== ##

getpoints <- function(y,estb0,sdb0=NULL,crit=qnorm(0.975))
{	
points(estb0,y,col="blue",pch=16)
if (!is.null(sdb0))
{
points(c(estb0-crit*sdb0,estb0+crit*sdb0),rep(y,2),col="red",type="l")
}
}
ordermethod <- c("gamma.ar1","logn.ar1","gamma.ind","logn.ind","semi.ar1","semi.ind")

estb0s <- c(
re.gamma.ar1$est[4,1],
re.logn.ar1$est[4,1],
re.gamma.ind$est[3,1],
re.logn.ind$est[3,1],
re.semi.ar1$est[4],
re.semi.ind$est[3]
)

## The true beta0 is:
b0
c <- 1.1
plot(0,0,type="n",xlim=c(min(estb0s)-0.5,max(estb0s)*c),ylim=c(0,7),yaxt="n",
main="Simulated from the AR(1) model \n with random effect ~ gamma")

legend("topright",
	legend=ordermethod)
abline(v=b0,lty=3)

## 1) gamma.ar1
sdb0 <- re.gamma.ar1$est[4,2]
getpoints(6,estb0s[1],sdb0)

## 2)logn.ar1
sdb0 <- re.logn.ar1$est[4,2] 
getpoints(5,estb0s[2],sdb0)

## 3) gamma.ind
sdb0 <- re.gamma.ind$est[3,2]
getpoints(4,estb0s[3],sdb0)

## 4) logn.ind
sdb0 <- re.logn.ind$est[3,2]
getpoints(3,estb0s[4],sdb0)

## 5) semi.ar1
getpoints(2,estb0s[5])

## 6) semi.ind
getpoints(1,estb0s[6])


}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

