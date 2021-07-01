\name{tvp.ridge}
\alias{tvp.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time-Varying Parameters via Ridge Regressions
%%  ~~function to do ... ~~
}
\description{This prototype function runs the 2-step Ridge Regression technique in Goulet Coulombe (2020). Basically, the technique allowing parameters of a linear equation(s) to evolve as driftless random walks and computations are carried via a closed-form solution in the form of data-augmented ridge regression rather than with typically burdensome filtering procedures. It handles the univariate and multivariate cases, and performs cross-validation of lambda (the amount of time variation). It will especially fast when T is below 1000, which is customary for macroeconomic time series, on which TVP models are often applied.
}


\usage{tvp.ridge(X,Y,lambda.candidates=exp(linspace(-6,20,n=15)),oosX=c(),
                      lambda2=.1,kfold=5,CV.plot=TRUE,CV.2SRR=TRUE,
                      sig.u.param=.75,sig.eps.param=.75,ols.prior=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Predictors (T * K)}
  \item{Y}{Targets (T * M), where M can be 1 (univariate case)}
  \item{lambda.candidates}{Candidate values for lambda (amount of time variation) to evaluate by CV.}
  \item{oosX}{A vector of X for which to produce one out-of-sample forecast. Since parameters evolve according to random walks, one should need to reevaluate at each t.}
  \item{lambda2}{Penalty on beta_0 (starting values of the random walks). Keep it soft unless X is large.}
\item{kfold}{Number of folds for CV (5 or 10, usually).}
\item{CV.plot}{If Switched to TRUE, the function will plot CV graph, which can be a useful check of whether your lambda.candidates are legit.}
\item{CV.2SRR}{Should we cross-validate again after the weight updating step? Can help marginally at low cost for small models. Can help marginally but incurs a higher computational cost for large X or Y.}
\item{sig.u.param}{How much should we shrink sigma_u estimates to the homogenous case? 1 is not at all, and 0 is completely. Intermediate values like 0.25, 0.5, and 0.75 are recommended. The formula is final.sig.u = [(hetero)^(sig.u.param)]*[(homo)^(1-sig.u.param)].}
\item{sig.u.param}{How much should we shrink sigma_eps estimates to the no-evolving volatility case? 1 is not at all, and 0 is completely. Intermediate values like 0.25, 0.5, and 0.75 are recommended. The formula is final.sig.eps = [(Garch(1,1))^(sig.u.param)]*[(homo)^(1-sig.u.param)].}
\item{ols.prior}{Will shrink beta_0 (at a strength given by lambda2) to OLS rather than 0. Not recommended unless X is very small. The Half & Half strategy described in the paper is typically preferable.}

%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{betas.rr}{Betas from a plain ridge regression (no two steps).}
\item{betas.2srr}{Betas from 2SRR (M * K * T).}
\item{lambdas}{Choice of lambda by CV. Will be a vector if M>1.}
\item{forecast}{A forecast for y if oosX was supplied. Will be a vector if M>1.}
\item{sig.eps}{Normalized evolving volatility weights.}
\item{yhat.rr}{In-sample span from a one-step approach. Will be a matrix if M>1.}
\item{yhat.2srr}{In-sample span from 2srr. Will be a matrix if M>1.}
}




\references{
Original paper is https://arxiv.org/abs/2009.00401}
\author{
Philippe Goulet Coulombe
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#data
data=matrix(rnorm(15*200),200,15)
X=data[,1:3]
Y=crossprod(t(X),rep(1,3))*(1-0.5*I(c(1:200)>75))+rnorm(200)/2

#2SRR
results <- tvp.ridge(Y=Y,X=X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
