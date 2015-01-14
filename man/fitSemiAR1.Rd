\name{fitSemiAR1}
\alias{fitSemiAR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the semi-parametric negative binomial mixed-effect AR(1) model.
}

\description{

%
This function fits the semi-parametric negative binomial mixed-effect AR(1) model in the formulation described Zhao et al (2013).
The conditional distribution of response counts given random effect is modelled by Negative Binomial as described in description of \code{\link{lmeNB}}. 
The conditional dependence among the response counts of a subject is modeled with AR(1) structure. 
The semiparametric procedure is employed for random effects. 
See descriptions of \code{\link{lmeNB}}.

}

\usage{

fitSemiAR1(formula, data, ID, Vcode,p.ini = NULL, IPRT = TRUE, deps = 0.001, maxit=100)

}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{formula}{
See \code{\link{lmeNB}}.
%
}
\item{data}{
See \code{\link{lmeNB}}.
}
%

\item{ID}{
See \code{\link{lmeNB}}.
}

\item{Vcode}{
See \code{\link{lmeNB}}.
}

  \item{p.ini}{
See \code{\link{fitParaAR1}}.
}
  \item{IPRT}{
See \code{\link{lmeNB}}.
}

  \item{deps}{
See \code{\link{lmeNB}}.
}

\item{maxit}{
See \code{\link{lmeNB}}.
}


}

\details{
The algorithm repeats the following four steps until a stoping criterion is satisfied:

Step 1) Estimate the coefficients of covariates by the method of weighted least squares.

Step 2) Approximate the distribution of the random effect \eqn{G_i}{G[i]} by \eqn{\gamma_i}{\gamma[i]}.

Step 3) Estimate \eqn{\alpha}{\alpha} and \eqn{\delta}{\delta} using the psudo-profile likelihood. 
This step calls \code{optim} to minimize the negative psudo log-likelihood with respect to \eqn{\log(\alpha)}{\log(\alpha)}) and logit(\eqn{\delta}{\delta}). The numerical integration is carried out using adaptive quadrature. When missing visits are present, the likelihood is approximated (See Zhao et al. 2013 for details).

Step 4) Estimate \eqn{Var(G_i)}{Var(G[i])} by the medhod of moment and update the weights.

All the computations are done in \code{R}.
}




%% \value{

%% \item{opt}{
%% The values returned by \code{optim} to minimize the negative of the psudo-profile log-likelihood with respect to \code{log(a)} and \code{logit(d)}.
%% }

%% \item{V}{
%% If the number of covariate is nonzero, \code{vcm} returns an naive estimate of the variance covariance matrix of \code{b0}, \code{b1},...
%% }

%% \item{est}{
%% A vector of length # covariates + 4 containing the estimates of the parameters \emph{(log(a), log(var(G)), logit(d), b0, b1, ...)}
%% }

%% %\item{gtb}{
%% %The relative frequency table of \emph{g_i}, (\emph{i=1,...,N}).
%% %\code{gh1} contains the unique values in ascending order and \code{ghw} contains relative frequencies.
%% %}

%% \item{gi}{
%% A vector of length \emph{N}, containing the approximation of \emph{g_i}s
%% }


%% \item{RE}{
%% \code{"NoN"}, denoting that the fitted model is a semi-parametric mixed-effect model.
%% }


%% %add
%% \item{Vcode}{A vector containing the input indices of time point.}

%% }

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
\code{\link{fitParaAR1}},
\code{\link{fitSemiIND}},
%\code{\link{fitSemiAR1}},

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
## ================================================================================ ##
## generate a data based on the semi-parametric negative binomial 
## mixed-effect AR(1) model.
## Under this model, the response counts follows the negative binomial:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## G_i is from unknown distribution.
## For simulation purpose, we generate the sample of gi from 
## the mixture of three gamma distribuions.

## The adjacent repeated measures of the same subjects are correlated 
## with correlation structure:
## cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)  

# log(a) = -0.5, log(th)=1.3, logit(delta) = -0.2
# b0 =  0.5, no covariates; 
loga <- -0.5
logtheta <- 1.3
logitd <- -0.2
b0 <- 0.5
# 80 subjects each with 5 scans
n <- 80
sn <- 5

## generate a sample of size B from the mixture of three gamma distribution:
p1 <- 0.5  
p2 <- 0.3
B <- 1000
sampledG<- c(
rgamma(n=p1*B,scale=1,shape=10),
rgamma(n=p2*B,scale=3,shape=5),
rgamma(n=(1-p1-p2)*B,scale=5,shape=5)
)


## mean is set to 1;
sampledG <- sampledG/mean(sampledG) 
logvarG <- log(var(sampledG))
## hist(sampledG)

DT4 <-  rNBME.R(gdist = "NoN",
               n = n, ## 	the total number of subjectss	       
	       sn = sn,
               u1 = rep(exp(b0),sn),
	       u2 = rep(exp(b0),sn),
	       a = exp(loga),
	       d = exp(logitd)/(1+exp(logitd)),
	       othrp = sampledG
	      )
Vcode<-rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
ID <- DT4$id
new <- Vcode > 0
dt4<-data.frame(CEL=DT4$y)
## ================================================================================ ##

## [1] Fit the negative binomial mixed-effect AR(1) model 
## where random effects is from the gamma distribution


re.gamma.ar1 <- fitParaAR1(formula=CEL~1,data=dt4,ID=ID,
		         Vcode=Vcode, 
		          p.ini=c(loga,logtheta,logitd,b0), 
		          ## log(a), log(theta), logit(d), b0
		          RE="G", 
		          IPRT=TRUE)

Psum<-index.batch(olmeNB=re.gamma.ar1, data=dt4,ID=ID,Vcode=Vcode,
	          labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE) 



## [2] Fit the negative binomial mixed-effect AR(1) model 
## where random effects is from the log-normal distribution


re.logn.ar1<-fitParaAR1(formula=CEL~1,data=dt4,ID=ID,
		        Vcode=Vcode, 
		        p.ini=c(loga,logtheta,logitd,b0), 
		        ## log(a), log(theta), logit(d), b0
		        RE="N", IPRT=TRUE)

## Requires some time
Psum<-index.batch(olmeNB=re.logn.ar1,data=dt4,ID=ID,Vcode=Vcode,
	          labelnp=new,qfun="sum", IPRT=TRUE) 



## [3] Fit the negative binomial independent model 
## where random effects is from the lognormal distribution
re.logn.ind<-fitParaIND(formula=CEL~1,data=dt4,ID=ID, 
                        RE="N", 			   	
		        p.ini=c(loga,logtheta,b0), 		
		        IPRT=TRUE)

Psum <- index.batch(olmeNB=re.logn.ind,data=dt4,ID=ID,
                    labelnp=new,qfun="sum", IPRT=TRUE) 


## [4] Fit the semi-parametric negative binomial AR(1) model 
## This model is closest to the true model

logvarG <- log(var(sampledG))

re.semi.ar1 <- fitSemiAR1(formula=CEL~1,data=dt4,ID=ID, 
                          p.ini=c(loga, logvarG, logitd,b0),Vcode=Vcode)
 
## compute the estimates of the conditional probabilities 
## with sum of the new repeated measure as a summary statistics 
Psum <- index.batch(olmeNB=re.semi.ar1, labelnp=new,data=dt4,ID=ID,Vcode=Vcode,
                    qfun="sum", IPRT=TRUE,i.se=TRUE) 

## compute the estimates of the conditional probabilities 
## with max of the new repeated measure as a summary statistics 
Pmax <- index.batch(olmeNB=re.semi.ar1, labelnp=new,qfun="max",data=dt4,ID=ID,Vcode=Vcode,
                    IPRT=TRUE,i.se=TRUE) 

## Which patient's estimated probabilities 
## based on the sum and max statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt4$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,1]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,1]


## [5] Fit the semi-parametric negative binomial independent model 


re.semi.ind <- fitSemiIND(formula=CEL~1,data=dt4,ID=ID, p.ini=c(loga, logvarG, b0))
Psum <- index.batch(olmeNB=re.semi.ind, labelnp=new,
                    data=dt4,ID=ID, qfun="sum", IPRT=TRUE,i.se=TRUE) 



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
ordermethod <- c("gamma.ar1","logn.ar1","logn.ind","semi.ar1","semi.ind")

estb0s <- c(
re.gamma.ar1$est[4,1],
re.logn.ar1$est[4,1],
re.logn.ind$est[3,1],
re.semi.ar1$est[4],
re.semi.ind$est[3]
)

## The true beta0 is:
b0
c <- 1.1
plot(0,0,type="n",xlim=c(min(estb0s)-0.5,max(estb0s)*c),ylim=c(0,7),yaxt="n",
main <- "Simulated from the AR(1) model \n with random effect ~ a semi-parametric distribution")

legend("topright",
	legend=ordermethod)
abline(v=b0,lty=3)

## [1] gamma.ar1
sdb0 <- re.gamma.ar1$est[4,2]
getpoints(6,estb0s[1],sdb0)

## [2] logn.ar1
sdb0 <- re.logn.ar1$est[4,2]
getpoints(5,estb0s[2],sdb0)

## [3] logn.ind
sdb0 <- re.logn.ind$est[3,2]
getpoints(4,estb0s[3],sdb0)

## [4] semi.ar1
getpoints(3,estb0s[4])

## [5] semi.ind
getpoints(2,estb0s[5])

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
