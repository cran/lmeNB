\name{index.batch}
\alias{index.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function to compute the point estimates and 95\% confidence intervals (for a parametric model) 
of the conditional probabilities 
\eqn{ 
Pr(q(\boldsymbol{Y}_{i,new}) \ge q(\boldsymbol{y}_{i,new})| \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{
Pr(q(Y[i,new])>=q(y[i,new])| Y[i,pre]=y[i,pre])
} for multiple subjects.
}

\description{
%
Let \eqn{m_i}{m[i]} be the number of pre-measurements and
\eqn{n_i}{n[i]} be the total number of repeated measures. 
Then the repeated measure of a subject can be divided into a pre-measurement set and a new measurement set as 
\eqn{
\boldsymbol{Y}_{i}=(\boldsymbol{Y}_{i,pre},\boldsymbol{Y}_{i,new})
}{
Y[i]=(Y[i,pre],Y[i,new])
}, where
\eqn{
\boldsymbol{Y}_{i,pre}=(Y_{i,1},\cdots,Y_{i,m_i})
}{
Y[i,pre]=(y[i,1],\cdots,Y[i,m[i]])
} 
and 
\eqn{
\boldsymbol{Y}_{i,new}=(Y_{i,m_i+1},\cdots,Y_{i,n_i})
}{
Y[i,new]=(Y[i,m[i]+1],...,Y[i,n[i]])
}
.
Given an output of \code{\link{fitParaIND}}, \code{\link{fitParaAR1}},
\code{\link{fitSemiIND}}, \code{\link{fitSemiAR1}} or \code{\link{lmeNB}},
this function computes the probability of observing the response counts as large as those new observations of subject \eqn{i}{i},
\eqn{
\boldsymbol{y}_{i,new}
}{
y[i,new]
} 
conditional on the subject's previous observations 
\eqn{
\boldsymbol{y}_{i,pre}
}{
y[i,pre]
} for subject \eqn{i}{i}.
That is, this function returns a point estimate and its asymptotic 95\% confidence interval (for a parametric model) of the conditional probability for each subject:
%

\eqn{ 
Pr(q(\boldsymbol{Y}_{i,new}) \ge q(\boldsymbol{y}_{i,new})| \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{ 
Pr(q(Y[i,new]) \ge q(y[i,new])| Y[i,pre]=y[i,pre]) 
}.

When the semiparametric approach is employed, the standard error and 95\% confidence intervals are computed using bootstrap samples.
A scalar statistic to summarize the new response counts can be either the total count, 
\eqn{q(\boldsymbol{Y}_{i,new})=\sum_{j=m_i+1}^{n_i} Y_{ij} }{  q(Y[i,new])=\sum[j=m[i]+1]^{n[i]} Y[ij]  }, 
or the maximum, 
\eqn{
q(\boldsymbol{Y}_{i,new})=\max\{ Y_{ij};j=m_i+1,\cdots,n_i \} 
}{
q(Y[i,new])=\max\{ Y[ij];j=m[i]+1,...,n[i] \} 
}.
%
See Zhao et al.(2013), for more details.
}



\usage{
index.batch(data,         labelnp,     ID,  Vcode = NULL, olmeNB = NULL, subset = NULL, 
            qfun = "sum", IPRT = TRUE, i.se = TRUE, MC = FALSE, C = FALSE,i.tol=1E-75)	
}


\arguments{
\item{data}{
See \code{\link{lmeNB}}.
This dataset does not have to be the same as the one used in the computations of negative binomial mixed effect regression 
(\code{\link{fitParaIND}}, \code{\link{fitParaAR1}},
\code{\link{fitSemiIND}}, \code{\link{fitSemiAR1}} or \code{\link{lmeNB}}).
}
\item{labelnp}{
A vector of length the total number of repeated measures 
(=\eqn{ \sum_{i=1}^N n_i }{ \sum[i=1]^N n[i] }), indicating new measures by \code{TRUE} and pre-measures by \code{FALSE}. For examples, suppose there are three subjects of interest.

The first subject has a \eqn{n_1=7}{n[1]=7} repeated measures and the last \emph{3} measures are new. 
The second and the third subjects both have \eqn{n_2=n_3=5}{n[2]=n[3]=5} repeated measures  and the last 2 measures are new. In this scenario,

\code{labelnp}\code{=c(rep(FALSE,4),} \code{rep(TRUE,3),} \code{rep(FALSE,3),} \code{rep(TRUE,2)} \code{,rep(FALSE,3)} \code{,rep(TRUE,2))}.
}
  \item{ID}{
   See \code{\link{lmeNB}}. The length of \code{ID} must be the same as \code{nrow(data)}.
 }
 \item{Vcode}{
    Necessary only if the \code{olmeNB} is an output of AR(1) models. 
    See \code{\link{lmeNB}}.
}   
  \item{olmeNB}{
    Output of \code{\link{fitParaIND}},\code{\link{fitParaAR1}},\code{\link{fitSemiIND}},\code{\link{fitSemiAR1}}
    or \code{\link{lmeNB}}.
  }
 \item{subset}{
   An optional expression indicating the subset of the subjects of that the index should be computed. 
 }
 
 \item{qfun}{
   If \code{qfun="sum"},
   a scalar statistic to summarize the new response counts is the total count.
   
   If \code{qfun="max"},
   a scalar statistic to summarize the new response counts is the maximum.
  
  %% linked by jCP_ar1
 }
 \item{IPRT}{
   print control.
 }
 
 
 \item{i.se}{
   If \code{i.se=TRUE} then the standard errors of the estimator of the conditional probability are returned for the output of 
   \code{\link{fitParaIND}} or \code{\link{fitParaAR1}}.
   The semi-parametric approach, \code{\link{fitSemiIND}} or \code{\link{fitSemiAR1}}, do not return the standard errors.
 }
 \item{MC}{
   Necessary when \code{olmeNB} if the AR(1) model outputs.%,\code{\link{fitParaAR1}} and \code{\link{fitSemiAR1}}.
   See \code{\link{CP.ar1.se}}.
 }

\item{C}{
See \code{\link{lmeNB}}. 
\code{C=TRUE} option could make computations of CPI faster for some patients.
}
\item{i.tol}{ See \code{\link{lmeNB}}.}
}


\details{
%
The standard error of the point estimate on the logit scale is constructed using the delta method for the parametric model, where distributional assumption was made for random effects. %put this in Details
%
}
\value{
The \eqn{N}{N} by 4 (3, if hide the SE) numeric matrix, containing the point estimate of the conditional probability, and the lower and the upper bounds of the 95 % confidence intervals. 
% hide the SE?
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

The internal functions of \code{\link{lmeNB}} for fitting relevant models:
\code{\link{fitParaIND}},
\code{\link{fitParaAR1}},
\code{\link{fitSemiIND}},
\code{\link{fitSemiAR1}},

The subroutines of \code{\link{index.batch}}:
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
## See the examples in help files of 
## fitParaIND, fitAR1IND, fitSemiIND, fitSemiAR1 and rNBME.R
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
