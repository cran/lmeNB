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
Zhao, Y., Li, D.K.B., Petkau, J.A., Riddehough, A. & Traboulsee, A.  Detection of unusual increases in MRI lesion counts in multiple sclerosis patients.
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
## generate a simulated dataset from the negative binomial mixed-effect 
## independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ log-normal(E(G_i)=1,var(G_i)=th)
set.seed(1)
sn = 5 # the number of repeated measures of each patient
n = 80 ## the number of patients
loga = - 0.5 
a = exp(loga) ## dispersion parameter 
logtheta <- 1.3
th = exp(logtheta) ## the variance of the gamma distributed random effect g_i


## No difference between the means of groups 
## The model only has an intercept term beta0 = 0.5
b0 = 0.5
u1 = rep(exp(b0),sn)  ## the mean vector for group 1 at time point 1,...,sn
u2 = rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn

## Data 0 generated from the IND model
DT.ind= rNBME.R(gdist="N", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn)
## Data 1
DT.ar= rNBME.R(gdist="N", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn, d=0.4)

dt.ind=data.frame(CEL=DT.ind$y,Vcode=DT.ind$vn-2,ID=DT.ind$id)
dt.ar=data.frame(CEL=DT.ar$y,Vcode=DT.ar$vn-2,ID=DT.ar$id)
## ================================================================================ ##

#### Fit IND models
tst1=lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=T)
tst2=lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=T)
tst3=lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=T, RE="N")
tst4=lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=T, RE="N") #not printing
tst5=lmeNB(CEL~Vcode, data=dt.ind, ID=dt.ind$ID, IPRT=T, RE="N")
## conditional probability index
Psum5=index.batch(olmeNB=tst5, labelnp=dt.ind$Vcode >= 1, data=dt.ind, ID=dt.ind$ID)



## Fit the semi-parametric model
tst6=lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=T, RE="semipara")
tst7=lmeNB(CEL~Vcode, data=dt.ind, ID=dt.ind$ID, IPRT=T, RE="semipara")
## conditional probability index
Psum7=index.batch(olmeNB=tst7,labelnp=dt.ind$Vcode >= 1, data=dt.ind, 
                 ID=dt.ind$ID, Vcode=Vcode)


## Fit the AR1 models
tst8=lmeNB(CEL~1, data=dt.ind, ID=dt.ind$ID, IPRT=T,AR=T,Vcode=dt.ind$Vcode)
tst9=lmeNB(CEL~1, data=dt.ar, ID=dt.ar$ID, IPRT=T,AR=T,Vcode=dt.ar$Vcode)
## conditional probability index
Psum9=index.batch(olmeNB=tst9, labelnp=dt.ind$Vcode >= 1, data=dt.ind, 
                  ID=dt.ind$ID,Vcode=dt.ind$Vcode)



}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
