
\name{lmeNB}
\alias{lmeNB}

%- Also NEED an '\alias' for EACH other topic documented here.
%
\title{
%
Performs the maximum likelihood estimation for the negative binomial
mixed-effect model.
This function is a wrapper for \code{\link{fitParaIND}},
\code{\link{fitParaAR1}}, \code{\link{fitSemiIND}} and \code{\link{fitSemiAR1}}.
}

\description{
Let \eqn{Y_{ij}}{Y[ij]} be the response count at \eqn{j}th repeated measure from \eqn{i}th subject.
The negative binomial mixed-effect independent model assumes that given the random effect \eqn{ G_i=g_i}{G[i]=g[i]}, 
the count response \eqn{Y_{ij}}{Y[ij]} follows the negative binomial distribution:
%

                   \eqn{ Y_{ij} | G_i =g_i \sim NB(r_{ij},p_i),}{Y[ij] | G[i] =g[i] ~ NB(r[ij],p[i]),} 

%
where  \eqn{p_i}{p[i]}, the failure probability of subject \eqn{i}{i} at each time point \eqn{j}{j} is parametrized as \eqn{p_i = \frac{1}{g_i\alpha+1} }{p[i]= 1/(g[i]*\alpha+1)} and \eqn{ \alpha > 0 }{\alpha>0}. The model assumes \eqn{ E(G_i) = 1 }{ E(G[i]) = 1} so that \eqn{ E(Y_{ij}|G_{i}=g_{i})=r_{ij}g_{i}\alpha }{E(Y[ij]|G[i]=g[i])=r[ij]*g[i]*\alpha} and \eqn{E(Y_{ij})= r_{ij}g_i}{E(Y[ij])= r[ij]*g[i]}. 
 This assumption allows the interpretation of the latent random variable \eqn{ G_i }{G[i]} as the subject i's activity rate relative to the overall cohort. The marginal mean \eqn{ \mu_{ij} = E(Y_{ij})}{ \mu[ij] = E(Y[ij])} is modeled with fixed effect coefficients, \eqn{\boldsymbol{\beta}}{\beta}: 
\eqn{\mu_{ij}= exp(\boldsymbol{X}_{ij}^T \boldsymbol{\beta})}{\mu[ij] = exp(X[ij]^T \beta)}.
Furthermore, let \eqn{Var(G_i)=\theta}{Var(G[i])=theta}, then \eqn{Var(Y_{ij})=\mu_{ij}^2\theta+\mu_{ij}(1+(\theta + 1)\alpha)}{Var(Y[ij])=mu[ij]^2theta+mu[ij](1+(theta + 1)alpha)}.


------------------------------------------------------------------------------------------------------------------

Regarding the dependence structures of \eqn{Y_{ij}}{Y[ij]} and \eqn{Y_{ij'}}{Y[ij']} conditional on the random effect \eqn{G_i}{G[i]}, we consider two models, namely independent and AR(1) models.
  
[1]: Independent model 
%
\eqn{Y_{ij}}{Y[ij]} and \eqn{Y_{ij'}}{Y[ij']} are independent conditionally on \eqn{G_i}{G[i]}. This assumption leads to 
%  		    
\eqn{
Cov(Y_{ij},Y_{ij'}|G_i=g_i)=0
}{
Cov(Y[ij],Y[ij']|G[i]=g[i])=0
} 
and 
\eqn{
Cov(Y_{ij},Y_{ij'})=\mu_{ij} \mu_{ij'}\theta
}{ 
Cov(Y[ij],Y[ij'])=mu[ij] mu[ij']theta
}
  
[2]: AR(1) model
%
Autoregressive (1) structures for \eqn{Y_{ij}}{Y[ij]} and \eqn{Y_ij'}{Y[ij']} conditionally on \eqn{G_i}{G[i]}.
That is given the random effect \eqn{ G_i=g_i}{G[i]=g[i]}, \eqn{Y_{ij}}{Y[ij]} depends on \eqn{Y_{i(j-1)}}{Y[i(j-1)]} 
through the beta binomial thinning and is conditionally independent on \eqn{Y_{ij'}}{Y[ij']} given 
\eqn{Y_{i(j-1)}}{Y[i(j-1)]} for all \eqn{j' < j-1 }{j' < j-1 }. 
The beta binomial thinning operator depends on a parameter \eqn{\delta}{\delta} and models to have 
\eqn{
Cov(Y_{ij},Y_{ij'}|G_i=g_i)=\delta^{j-j'} \mu_{ij'}(\alpha g_i^2+g_i)
}{
Cov(Y[ij],Y[ij']|G[i]=g[i])=\delta^{j-j'} \mu[ij](\alpha*g[i]^2+g[i])
}.
This means that \eqn{\delta}{\delta} measures the strength of the positive AR(1) association between repeated measures of a subject: 
the larger \eqn{\delta}{\delta}, the stronger the positive correlations between the repeated measures of the same subject are.

-----------------------------------------------------------------------------------------------------------------
  
Regarding the random effect \eqn{G_i}{G[i]} distribution, \code{lmeNB} allows three models, namely log-normal, gamma and semiparametric models.
(All models assume  \eqn{E(G_i)=1}{E(G[i])=1} and \eqn{Var(G_i)=\theta}{Var(G[i])=\theta}.)

(1) The log-normal model
%
That is regular log-normal parameters are restricted as 
\code{meanlog=-log(theta+1)/2}, 
\code{sdlog = sqrt(log(theta+1))}.

(2) The gamma model
%
That is regular gamma parameters are restricted as \code{shape=1/theta}, \code{scale=theta}.
  
(3) The semiparametric model
%
No distributional assumption and the random effect distribution is approximated by the estimated values of the quantity:
\eqn{
\gamma_i = w_i \frac{ y_{i+} }{ \mu_{i+}} + (1-w_i)
}{
\gamma[i] = w[i] (y[i+]/ \mu[i+]) + (1 - w[i])
},
where:
\eqn{
y_{i+} = \sum_{j=1}^{n_i} y_{ij} 
}{
y[i+] = \sum[j=1]^{n[i]} y[ij] 
}, 
\eqn{
\mu_{i+} =  \sum_{j=1}^{n_i} \mu_{ij} 
}{
\mu[i+] =  \sum[j=1]^{n[i]} \mu[ij] 
} 
and, 
\eqn{
w_{i} = \sqrt{ \frac{Var(G_i) }{ Var(Y_{i+}/\mu_{i+}) } }
}{
w[i] = sqrt( Var(G[i])/Var(Y[i+]/\mu[i+]) )
}.
See Zhao et al. (2013) for more details.
%
The missing count responses, if assumed to be missing at random, can be igored. Other types of missing data are currently not accepted. 
  %  


}


\usage{
lmeNB(formula,data,ID,p.ini=NULL,IPRT=FALSE,AR=FALSE,
      RE=c("G","N","NoN"),deps=1e-03,Vcode=NULL,C=FALSE,
      i.tol=1e-7,o.tol=sqrt(.Machine$double.eps),labelnp,
      maxit=100,semi.boot=0,u.low=0)
}

\arguments{

\item{formula}{
An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
The formula must contain an intercept term.
}
\item{data}{
A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
The each row must contains the data corresponding to the repeated measure \eqn{j}{j} of subjects and the rows \eqn{(i,j)}{(i,j)}s must be ordered in a way that measurements from a subject is clustered together as \eqn{(1,1)}{(1,1)},...,\eqn{(1,n_1)}{(1,n[1])},\eqn{(2,1)}{(2,1)},...,\eqn{(2,n_2)}{(2,n[2])},...,\eqn{(N,n_N)}{(N,n[N])}. Missing values are accepted of the response variables are treated as missing at random and simply removed from the data when \code{AR=FALSE}. When \code{AR=TRUE}, then the missing values are allowed. See the reference for missing value treatments. Missing values in covariates are not currently accepted.
}
\item{ID}{
A vector of length \eqn{ \sum_{i=1}^N n_i }{ \sum[i=1]^N n[i] }, containing the patient IDs that corresponds to \code{data}. i.e., \code{c(rep(ID_1,n_1),}\code{rep(ID_2,n_2)}\code{,...,rep(ID_N,n_N))}. The length must be the same as the number of rows of \code{data}. Missing ID values are NOT accepted.
}
  %
\item{p.ini}{
The initial values of the parameters:

If \code{AR=0},
\code{p.ini}=
\eqn{
(\log(\alpha), \log(Var(G_i))
}{
(\log(\alpha), \log(Var(G[i]))}
\eqn{
,\beta_0, \beta_1,...)
}{
,beta[0], beta[1],...)
}.

If \code{AR=1},
\code{p.ini}=\eqn{
(\log(\alpha), \log(Var(G_i)),logit(\delta),\beta_0,\beta_1,...)
}{
(\log(\alpha), \log(Var(G[i])),logit(\delta),\beta[0], \beta[1],...)
}.

\code{NULL} is accepted.
}

\item{IPRT}{
A logical, passed to Iprt of function \code{optim}. If \code{TRUE} then print iterations. 
}

\item{AR}{
  A logical, if \code{TRUE}, then the AR(1) model is employed. 
  If \code{FALSE}, then the independent model is employed.
}
  
\item{RE}{
  The distribution of random effects \eqn{G_i}{G[i]}.
  If \code{RE="G"}, then the conditional probability is computed by assuming the random effect is from a gamma distribution with mean 1 and variance \eqn{\theta}{theta} (gamma model).
  If \code{RE="N"}, then the conditional probability is computed by assuming the random effect is from a log-normal distribution with mean 1 and variance \eqn{\theta}{theta} (log-normal model).
  If \code{RE="NoN"}, then the conditional probability is computed based on the semi-parametric model with mean 1 (semiparametric model).
  % refered by jCP_ar1
}
\item{deps}{
 In the semiparametric models, the algorithms are terminated when the maximum difference of fixed effect coefficients between the current and previous iterations is less than \code{deps}. Passed to \code{\link{fitSemiIND}} and \code{\link{fitSemiAR1}}.
}

\item{Vcode}{
Necessary only if the AR(1) model is fit.
A vector of length the total number of repeated measures, containing the indices of time point.
For example, there are three subjects and the first two subjects do not have missing visits and completed five visits while the last  subject missed the third visit  and have four visits in total, then
\code{Vcode=c(1,2,3,4,5,1,2,3,4,5,1,2,4,5)}. 
}






\item{i.tol}{
The absolute tolerance of \code{\link{integrate}} function, which is used to integrate out the random effect of every patient.
Used only in parametric methods.
\code{i.tol} should be about \code{1E-3} for \code{C=TRUE} option.
}
\item{o.tol}{
The relative tolerance for \code{\link{optim}} function which is used to search for the maximum likelihood.
Used only in parametric methods.
}
\item{labelnp}{
See \code{\link{index.batch}}.
Necessary only for semiparametric random effect model and \code{semi.boot > 1}.
To account for the varying follow-up times of the patients, the bootstrap sampling is stratified according to the follow-up time.
}
\item{maxit}{
  The maximum number of iterations.
  Necessary only for semiparametric random effect model.
%%Write here !! Cite SIM 
}

\item{C}{
If \code{C=TRUE}, then the function uses the likelihood written in \code{C}.
The integration of the random effect is done using Cubature (Multi-dimensional integration) package written by Steven G. Johnson.
This option could make computation faster in some scenario.
If \code{C=FALSE}, then the function the likelihood is likelihood written in R language.
The integration of the random effect is done using \code{integrate} function in the \code{stats} package.      
}
\item{semi.boot}{
The number of bootstrap samples to construct the bootstrap empirical confidence intervals for the fixed effect coefficients. If the value is less than 1 then the bootstrap confidence intervals are not computed. Necessary only for semiparametric random effect model.
}

\item{u.low}{
For the semiparametric procedures, 
we observed that the algorithm could behave very unstable when factor covariates are employed in the dataset, 
and data contains "few" repeated measures of one of the corresponding factor groups:
The algorithm could yield unacceptably small estimate of \eqn{\mu_{ij}}{mu[ij]}.% that machine treats as zero. 
As 
\eqn{Var(Y_{ij})}{Var(Y[ij])} and 
\eqn{Cov(Y_{ij},Y_{ij'})}{Cov(Y_{ij},Y_{ij'})} are both multiples of \eqn{\mu_{ij}}{mu[ij]} (see description above), 
small \eqn{\mu_{ij}}{mu[ij]} leads to singular \eqn{Var(\boldsymbol{Y}_{i})}{Var(Y[i])} matrix. 
As 
a result, the algorithm breakdown when computing the weighted matrix of the weighted least square,
\eqn{\boldsymbol{W}_i}{W[i]}=\eqn{Var(\boldsymbol{Y}_{i})^{-1}}{Var(Y[i])^{-1}}.
To prevent this issue, the current algorithm takes add-hoc treatment, and replaces small \eqn{\mu_{ij}}{mu[ij]} i.e.
those with \eqn{\mu_{ij}<}{mu[ij]<} \code{u.low} with \code{u.low}
when calculating the weight matrix \eqn{\boldsymbol{W}_i}{W[i]}.
%
\code{u.low}=0 means that there is no modification, and the smaller \code{u.low}, the "closer" the modified algorithm is to the original one proposed in Zhao et al. (2013). Our preliminary study indicates that \code{u.low> 1E-4} prevents the breakdown problem and the performance of the algorithm are similar when \code{1E-3 < u.low < 1E-1} in terms of the root mean square error of the conditional probability index.  

}





}
\details{

\code{fitParaIND} calls \code{optim} to minimize the negative log-likelihood of 
the negative binomial model with respect to the model parameters: 
\eqn{
\log(\alpha), \log(\theta),\beta_0, \beta_1,...
}{
\log(\alpha), \log(\theta),\beta[0], \beta[1],...
}. 
The Nelder-Mead algorithm is employed.
The log-likelihood is obtained by marginalizing out the random effects. 
The numerical integration is carried out using adaptive quadrature.
The missing count responses, if assumed to be missing at random, can be igored.
Other types of missing data are currently not accepted.




When the estimated over-dispersion parameter (\eqn{\alpha}{\alpha}) is close to zero, the negative binomial model reduces to the poisson model, suggesting that the negative binomial mixed-effect model might not be appropriate.
When \code{AR}=1 and the estimated auto-correlation parameter (\eqn{\delta}{\delta}) is close to zero, the model is suggesting that there is no AR(1) structure among the sequentially collected responses. Hence user might use \code{AR}=0 setting which assumes no AR(1) structure.
%


We note that the results could be sensitive to initial values.
}

\value{



\item{call}{
The input of the function.
}

\item{opt}{

If \code{RE="G"} or \code{RE="N"},
then \code{opt} contains the results directly from the \code{optim} function, 
which is used to minimize the negative of the log-likelihood. 

If \code{RE="NoN"},
then \code{opt} contains the results directly from the \code{optimize} function,
which is used to minimize the negative of the psudo-profile log-likelihood with respect to the dispersion parameter alpha, \code{a}.
}

\item{nlk}{
The value of the negative log-likelihood corresponding to \code{opt$par}. 
}

\item{V}{ 
If \code{RE="G"} or \code{RE="N"},
the approximated asymptotic covariance matrix of the maximum likelihood estimators, i.e. \code{V}=\code{solve(opt$hessian)}.
If \code{opt$hessian} is not invertible, then \code{V} = \code{matrix(NA,length(p.ini),length(p.ini))}.

If \code{RE="NoN"}, \code{AR=FALSE} and \code{semi.boot>0},
\code{V} contains the bootstrap covariance matrix based on \code{semi.boot} samples. 
}

\item{est}{ 
The matrix of the number of fixed-effect parameters (i.e. \code{length(p.ini)}) by 2. 
%%
The first column contains the estimates of the model parameters.
%
The second column contains the standard error of the estimators, i.e., \code{sqrt(diag(V))}.
%
If \code{V} is not evaluated, then \code{est} only has one column.

}



}

\references{
Detection of unusual increases in MRI lesion counts in individual multiple sclerosis patients. (2013) Zhao, Y., Li, D.K.B., Petkau, A.J., Riddehough, A., Traboulsee, A., Journal of the American Statistical Association.

A flexible mixed effect negative binomial regression model for detecting unusual increases in MRI lesion counts in individual multiple sclerosis patients. Kondo, Y., Zhao, Y., Petkau, A.J.
}
\author{
 Zhao, Y. and Kondo, Y.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The subroutines of this function is:
\code{\link{fitParaIND}},
\code{\link{fitParaAR1}},
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
## ================================================================================ ##
## generate a simulated dataset from the negative binomial mixed-effect 
## independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ log-normal(E(G_i)=1,var(G_i)=th)
set.seed(1)
sn <- 5 # the number of repeated measures of each patient
n <- 80 ## the number of patients
loga <- - 0.5 
a <- exp(loga) ## dispersion parameter 
logtheta <- 1.3
th <- exp(logtheta) ## the variance of the gamma distributed random effect g_i


## No difference between the means of groups 
## The model only has an intercept term beta0 <- 0.5
b0 <- 0.5
u1 <- rep(exp(b0),sn)  ## the mean vector for group 1 at time point 1,...,sn
u2 <- rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn

## DT.ind is generated from the IND model
DT.ind<- rNBME.R(gdist="N", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn)
## DT.ar is generated from AR(1) model with delta=0.4
DT.ar<- rNBME.R(gdist="N", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn, d=0.4)

dt.ind<-data.frame(CEL=DT.ind$y,Vcode=DT.ind$vn-2,ID=DT.ind$id)
dt.ar<-data.frame(CEL=DT.ar$y,Vcode=DT.ar$vn-2,ID=DT.ar$id)
## ================================================================================ ##
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#### Fit various parametric independent models ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
tst1 <- lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE)            ## gamma Gi
tst2 <- lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=TRUE)              ## gamma Gi
tst3 <- lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE, RE="N")    ## log-normal Gi
tst4 <- lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=TRUE, RE="N")      ## log-normal Gi
tst5 <- lmeNB(CEL~Vcode, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE, RE="N")## log-normal Gi
## conditional probability index with the fitted results of tst5
Psum5 <- index.batch(olmeNB=tst5, labelnp=dt.ind$Vcode >= 1, data=dt.ind, ID=dt.ind$ID)


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#### Fit the semi-parametric independent model ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
tst6 <- lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE, RE="NoN")
tst7 <- lmeNB(CEL~Vcode, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE, RE="NoN",
              semi.boot=100,labelnp=dt.ind$Vcode >= 1) 
## semi.boot = 100 option computes bootstrap SE and 95%CI based on 100 bootstrap samples 
## conditional probability index with fitting result of tst7
Psum7 <- index.batch(olmeNB=tst7,labelnp=dt.ind$Vcode >= 1, data=dt.ind, 
                     ID=dt.ind$ID, Vcode=Vcode)

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#### Fit the parametric AR1 models ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
tst8 <- lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=TRUE,AR=TRUE,Vcode=dt.ind$Vcode)
tst9 <- lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=TRUE,AR=TRUE,Vcode=dt.ar$Vcode)
## conditional probability index
Psum9 <- index.batch(olmeNB=tst9, labelnp=dt.ind$Vcode >= 1, data=dt.ind, 
                     ID=dt.ind$ID,Vcode=dt.ind$Vcode)


}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
