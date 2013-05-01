\name{rNBME.R}
\alias{rNBME.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a dataset from the negative binomial mixed-effect independent/AR(1) model 
}
\description{
This function simulates a dataset based on the negative binomial mixed-effect independent/AR(1) model with two treatment groups described in Zhao et al. The group mean can be different at each time point, but no other covariates are allowed.
%
See \code{\link{mle.fun}}, \code{\link{mle.ar1.fun}} for details of the model explanations.
%
}
\usage{
rNBME.R(
	gdist = "G", n = 200, sn = 5, th = exp(1.3), 
	u1 = rep(1.5, 5), u2 = rep(1.5, 5), 
	a = exp(-0.5),d=NULL,  othrp = list(u.n = 3, s.n = 0.5, p.mx = 0.05, sh.mx = NA)
	)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gdist}{

The distribution of the random effect term \emph{G_i}.

If \code{gdist="G"}, \emph{G_i} is from the gamma distribution.

If \code{gdist="N"},  \emph{G_i} is from the log normal distribution.

If \code{gdist="U"}, \emph{G_i} (on the log scale) is from the uniform distribution.

If \code{gdist="GN"}, \emph{G_i} is from the mixture of the gamma distribution and the normal distribution.

If the generated values are negative, they are truncated to zero.

If \code{gdist="NoN"}, \emph{G_i} is sampled from the pre-specified vector \code{othrp} with replacement.




}
  \item{n}{
The number of patients. It must be an even number.
}
  \item{sn}{
The number of repeated measures per patient. 
Generated datasets are balanced design.    
}
  \item{th}{
If \code{gdist="G"}, \code{th} is a \code{scale} parameter of the gamma distribution.

If \code{gdist="N"} or \code{gdist=="U"}, \code{th} is \emph{var(G_i)}.

If \code{gdist="GN"}, see details.

If \code{gdist="NoN"}, this parameter is not used.
} 

\item{u1}{
%
A vector of length \emph{sn}, specifying the mean of the treatment group 1 \emph{ E(Y_ij)} =\code{ u1[j]}.
%
% Although \emph{beta} of the negative binomial mixed-effect independent model does not depend on time points,
%
%\code{ind.sdt} allows to generate response counts whose mean depend on time %points.
%
}
\item{u2}{
A vector of length \emph{sn}, specifying the mean of the treatment group 2 \emph{ E(Y_ij)} =\code{ u2[j]}.
}
  \item{a}{
%
The parameter \emph{a} of the negative binomial mixed-effect independent model. See \cite{\link{mle.fun}}. 
%
}

\item{d}{
If \code{d=NULL}, generate data from the independent model.
If \code{d} is a scalar between 0 and 1, then d is \emph{delta} in the AR(1) model, and generate datasets from the AR(1) model.
}
  \item{othrp}{
If \code{gdist="GN"}, parameters for the \code{GN} option. See details.
If \code{gdist="NoN"}, \code{othrp} is a vector, containing a sample of \emph{G_i}, which is treated as a population and \emph{G_i} is resampled.
}
}


%
%
\details{
%
The generated datasets have equal number of scans per person.
%
The number of patients in the two groups are the same.
%
If \code{gdist=="GN"}, datasets are generated from:
%


\code{othrp$p.mx}*N(\code{mean}=\code{othrp$u.n},\code{s.d}=\code{othrp$s.n}) + (1-\code{othrp$p.mx})*gamma\code{(scale=th,shape)},
%
where \code{shape} of the gamma distribution is chosen to ensure  \emph{E(G_i)=1}.


}

\value{

\item{id}{The vector of length \code{n*sn} containing patient IDs: \code{rep(1:n,each=sn)}}

\item{vn}{
The vector of length \code{n*sn} containing the indicies of time points: \code{rep(1:sn, n)}}
\item{gp}{The vector of length \code{n*sn} containing the indicies of the treatment groups}
\item{y}{The vector of length \code{n*sn} containing generated response counts} 
\item{g}{The vector of length \code{n*sn} containing generated random effect terms}
\item{Gpara}{The record of the distribution and parameter specifications used to generate the dataset} 

}

\references{
 Zhao, Y., Li, D.K.B., Petkau, J.A., Riddehough, A. Traboulsee, A.  Detection of unusually large increases in 
MRI lesion counts in multiple sclerosis patients.
}
\author{
 Zhao, Y. and Kondo, Y.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
The functions to fit related models:
\code{\link{mle.fun}},
\code{\link{mle.ar1.fun}},
\code{\link{mle.a3.fun}},
\code{\link{mle.ar1.non3}},

The subroutines of \code{\link{index.batch}} to compute the conditional probability index:
\code{\link{jCP.ar1}},
\code{\link{CP1.ar1}},
\code{\link{MCCP.ar1}},
\code{\link{CP.ar1.se}},
\code{\link{CP.se}},
\code{\link{jCP}},

%The functions to generate simulated datasets from the negative binomial mixed-effect AR(1) model:
%\code{\link{rNBME.R}}.
}
\examples{

\dontrun{

## ================================================================================ ##
## generate a simulated dataset from the negative binomial mixed-effect independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ log-normal(E(G_i)=1,var(G_i)=th)
set.seed(1)
sn = 5 ## the number of repeated measures of each patient
n = 80 ## the number of patients
loga = - 0.5
a = exp(loga) ## the parameter for the failure probability of the negative binomial distribution
logtheta <- 1.3
th = exp(logtheta) ## the parameter for the gamma distribution of g_i

## No difference between the means of groups 
## The model only has an intercept term beta0 = 0.5
b0 = 0.5
u1 = rep(exp(b0),sn) ## the mean vector for group 1 at time point 1,...,sn
u2 = rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn


DT0= rNBME.R(gdist="N", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn)
ID = DT0$id
Vcode=rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
new = Vcode > 0
dt5=data.frame(CEL=DT0$y)
logitd = -0.5

## ================================================================================ ##

## 1) Fit the negative binomial mixed-effect AR(1) model where random effects are from the gamma distribution


re.gamma.ar1=mle.ar1.fun(formula=CEL~1,data=dt5,ID=ID,
		      Vcode=Vcode, 
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="G", 
		      IPRT=TRUE)

Psum=index.batch(obj=re.gamma.ar1, 
	         ct=new,qfun="sum", iprt=TRUE,i.se=FALSE) ## i.se=T requires more time...
		 


## 2) Fit the negative binomial mixed-effect AR(1) model where random effects is from the log-normal distribution

re.logn.ar1=mle.ar1.fun(formula=CEL~1,data=dt5,ID=ID,
		      Vcode=Vcode, 
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="N", 
		      IPRT=TRUE)

Psum=index.batch(obj=re.logn.ar1,
	         ct=new,qfun="sum", iprt=TRUE) 



## 3) Fit the negative binomial independent model where random effects is from the gamma distribution
## This is the true model

re.gamma.ind=mle.fun(formula=CEL~1,data=dt5,ID=ID, 
                   model="G", 
		   p.ini=c(loga,logtheta,b0), 
		   IPRT=TRUE)

## compute the conditional probability with sum of the new repeated measure as a summary statistics 
Psum=index.batch(obj=re.gamma.ind, 
                 ct=new,qfun="sum", iprt=TRUE) 

## compute the conditional probability with max of the new repeated measure as a summary statistics 
Pmax=index.batch(obj=re.gamma.ind, 
                 ct=new,qfun="max", iprt=TRUE) 

## Which patient's estimated probabilities based on the sum and max statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax[,1]-Psum[,1]))==n) )
## Show the patient's CEL counts  
dt5$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax[IDBigDif,]



## 4) Fit the negative binomial independent model where random effects is from the lognormal distribution
re.logn.ind=mle.fun(formula=CEL~1,data=dt5,ID=ID, 
                   model="N", 			   	
		   p.ini=c(loga,logtheta,b0), 		
		   IPRT=TRUE)

Psum=index.batch(obj=re.logn.ind, 
                 ct=new,qfun="sum", iprt=TRUE) 


## 5) Fit the semi-parametric negative binomial AR(1) model 

logvarG = -0.5

re.semi.ar1=mle.ar1.non3(formula=CEL~1,data=dt5,ID=ID, Vcode=Vcode,p.ini=c(loga, logvarG, logitd,b0))
Psum=index.batch(obj=re.semi.ar1, ct=new,qfun="sum", iprt=TRUE)  



## 6) Fit the semi-parametric negative binomial independent model 
re.semi.ind=mle.a3.fun(formula=CEL~1,data=dt5,ID=ID, p.ini=c(loga, logvarG, b0))
Psum=index.batch(obj=re.semi.ind,ct=new, qfun="sum", iprt=TRUE) 



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
main="Simulated from the independent model \n with random effect ~ log-normalnbinomRE_1.0.zip")

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
