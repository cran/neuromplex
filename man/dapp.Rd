\name{dapp}
\Rdversion{0.5}
\alias{dapp}
\title{Dynamic Admixture of Poisson Process}
\description{Fits the DAPP model to binned spiking data}
\usage{ 
dapp(spike.counts, lengthScale = NULL, lsPrior = NULL,
     hyper = list(prec = c(1,1), sig0 = 1.87, w=c(1,1)),
     burnIn = 1e3, nsamp = 1e3, thin = 4,
     verbose = TRUE, remove.zeros = FALSE)
} 
\arguments{
\item{spike.counts}{A list with the following items. 'Acounts': binned spike counts under condition A presented as a matrix. Rows are bins, columns are replicates (trials). 'Bcount': binned spike counts under condition B. 'ABcounts': binned spike counts under condition AB. 'bin.mids': an array giving the mid-points of the time bins. 'bin.width': a scalar giving the bin width.}
\item{lengthScale}{an array giving the length scale parameter values to be used for Gaussian process prior. Defaults to \code{sort(0.16 * resp.horiz / c(4, 3, 2, 1, 0.5, 0.1))} where \code{resp.horiz} is the time horizon of the response period.}
\item{lsPrior}{an array of the same length as \code{lengthScale} giving the prior probabilities of the length scale values.}
\item{hyper}{a list of hyper parameters with the following iterms. 'prec': a 2-vector giving the shape and rate parameters of the gamma distribution on the Dirichlet precision parameter. 'sig0': a scalaer giving the scale of the (centered) logistic distribution used in transforming the Gaussian random curves into curves restricted between 0 and 1.}
\item{burnIn}{number of MCMC iterations to discard as burn-in.}
\item{nsamp}{number of MCMC draws to be saved for posterior inference.}
\item{thin}{the thinning rate at which MCMC draws are to be saved. The total number of iterations equals \code{burnIn + nsamp * thin}}
\item{verbose}{logical indicating if some fit details should be printed during the course of the MCMC}
\item{remove.zeros}{logical indicating if trials with zero spike count shuold be removed from the analysis}
}
\value{
Returns a list of class "dapp" containting the following items.

\item{lsProb}{posterior preditctive draws of length scale}
\item{lambda.A}{posterior draws of lambda.A at bin mid-points}
\item{lambda.B}{posterior draws of lambda.B at bin mid-points}
\item{alpha}{posterior draws of the alpha curves at bin mid-points}
\item{A}{posterior draws of the latent variable A which gives the AB spike counts (by bin) that are to be attributed to signal A (the remaining are attributed to signal B)}
\item{prec}{posterior draws of precision}
\item{alpha.pred}{posterior predictive draws of alpha (of a future trial)}
\item{psl.pred}{posterior predictive draw of the feature parameters (phi, psi, ell) (of a future trial)}
\item{details}{mcmc details given as an array of \code{c(niter, nsamp, burnIn, thin, MH acceptance rate)}}
\item{hyper}{hyper parameters used in model fitting}
\item{lengthScale}{length scale set used in model fitting}
\item{lsPrior}{length scale prior}
\item{bin.mids}{bin mid-points}
\item{bin.width}{bin width}
\item{mcmc}{mcmc controls (burn-in length, thinning rate and number of saved draws)}
}

\references{Glynn, C., Tokdar, S.T., Zaman, A., Caruso, V.C., Mohl, J.T., Willett, S.M., and Groh, J.M. (2020+). Analyzing second order stochasticity of neural spiking under stimuli-bundle exposure. The Annals of Applied Statistics. Accepted.}
\seealso{\code{\link{plot.dapp}}, \code{\link{summary.dapp}} and \code{\link{predict.dapp}}.}

\examples{
  ## Note:
  #### The example below uses a simpler synthetic data, a wider bin-width
  #### and a shorter MCMC run to keep the run length less than 5s
  #### Use ?plot.dapp or ?plot.summary for a more realistic example
  
  ## Generate 25 A and 30 B trials with rate functions
  ##    lambda.A(t) = 160*exp(-2*t/1000) + 40*exp(-0.2*t/1000)
  ##    lambda.B(t) = 40*exp(-2*t/1000)
  ## where time t is measured in ms. Then, generate 30 AB trials,
  ## roughly 2/3 with flat weight curves with a constant intensity
  ## either close to A, or close to B (equally likely). The 
  ## remaining 1/3 curves are sinusoidal that snake between 0.01 and 0.99 
  ## with a period randomly drawn between 400 and 1000
  
  ntrials <- c(nA=25, nB=30, nAB=30)
  flat.range <- list(A=c(0.85, 0.95),
                     B=c(0.05, 0.15))
  flat.mix <- c(A=1/2, B=1/2)
  wavy.span <- c(0.01, 0.99)
  wavy.period <- c(400, 1000)
  
  T.horiz <- 1000
  rateB <- 40 * exp(-2*(1:T.horiz)/T.horiz)
  rateA <- 4*rateB + 40 * exp(-0.2*(1:T.horiz)/T.horiz)
  
  synth.data <- synthesis.dapp(ntrials = ntrials, pr.flat = 2/3,
                               intervals = flat.range, wts = flat.mix,
                               span = wavy.span, period.range = wavy.period,
                               lambda.A=rateA, lambda.B=rateB)
  
  ## Generate binned spike counts witb 100 ms bins
  spike.counts <- mplex.preprocess(synth.data$spiketimes, bw=100, visualize=FALSE)
  
  ## Fit the DAPP model to data
  #### A short MCMC run is done below to keep the run length short.
  #### Use default or larger values for burn, nsamp and thin
  #### for more reliable estimation
  fit.post <- dapp(spike.counts, burn=10, nsamp=90, thin=1, verbose=FALSE)
  
  ## Visualize model fit
  plot(fit.post)
  
  ## Post process results to assign second order stochasticity labels
  summary(fit.post)
}

\keyword{programming}
