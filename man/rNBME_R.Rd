\name{rNBME.R}
\alias{rNBME.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a dataset from the negative binomial mixed-effect independent/AR(1) model 
}
\description{
This function simulates a dataset based on the negative binomial mixed-effect independent/AR(1) model with two treatment groups described in Zhao et al (2013). The group mean can be different at each time point, but no other covariates are allowed.
%
See \code{\link{fitParaIND}}, \code{\link{fitParaAR1}} for details of the model explanations.
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

The distribution of the random effect term \eqn{G_i}{G[i]}.

If \code{gdist="G"}, \eqn{G_i}{G[i]} is from the gamma distribution.

If \code{gdist="N"}, \eqn{G_i}{G[i]} is from the log normal distribution.

If \code{gdist="U"}, \eqn{G_i}{G[i]} (on the log scale) is from the uniform distribution.

If \code{gdist="GN"}, \eqn{G_i}{G[i]} is from the mixture of the gamma distribution and the normal distribution.

If the generated values are negative, they are truncated to zero.

If \code{gdist="NoN"}, \eqn{G_i}{G[i]} is sampled from the pre-specified vector \code{othrp} with replacement.

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

If \code{gdist="N"} or \code{gdist=="U"}, \code{th} is \eqn{Var(G_i)}{Var(G[i])}.

If \code{gdist="GN"}, see details.

If \code{gdist="NoN"}, this parameter is not used.
} 

\item{u1}{
%
A vector of length \emph{sn}, specifying the mean of the treatment group 1 \eqn{ E(Y_ij)}{E(Y[ij])} =\code{ u1[j]}.
% Although \emph{beta} of the negative binomial mixed-effect independent model does not depend on time points,
%\code{ind.sdt} allows to generate response counts whose mean depend on time %points.
%
}
\item{u2}{
A vector of length \eqn{sn}{sn}, specifying the mean of the treatment group 2 \eqn{ E(Y_ij)}{E(Y[ij])} =\code{ u2[j]}.
}
  \item{a}{
%
The dispersion parameter \eqn{\alpha}{\alpha} of the negative binomial mixed-effect independent model. 
See description in \cite{\link{lmeNB}}. 
%
}

\item{d}{
If \code{d=NULL}, generate data from the independent model.
If \code{d} is a scalar between 0 and 1, then d is \eqn{\delta}{delta} in the AR(1) model, and generate datasets from the AR(1) model.
}
  \item{othrp}{
If \code{gdist="GN"}, parameters for the \code{GN} option. See details.
If \code{gdist="NoN"}, \code{othrp} is a vector, containing a sample of \eqn{G_i}{G[i]}, which is treated as a population and \eqn{G_i}{G[i]} is resampled.
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
where \code{shape} of the gamma distribution is chosen to ensure  \eqn{E(G_i)=1}{E(G[i])=1}.


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

The functions to fit related models:
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

%The functions to generate simulated datasets from the negative binomial mixed-effect AR(1) model:
%\code{\link{rNBME.R}}.
}
\examples{

\dontrun{
## See the examples in help files of fitParaIND, fitParaAR1, fitSemiIND, fitSemiAR1 and lmeNB

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
