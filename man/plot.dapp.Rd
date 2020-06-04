\name{plot.dapp}
\Rdversion{0.5}
\alias{plot.dapp}
\title{Plotting Method for Dynamic Admixture of Poisson Process}
\description{Visually summarizes model fit of the DAPP model to binned spiking data}
\usage{ 
\method{plot}{dapp}(x, tilt.prior = FALSE, mesh.tilt = 0.1,
     nprior = x$mcmc["nsamp"], ncurves = 10,
     simple.layout = FALSE, ...)
} 
\arguments{
\item{x}{a fitted model of the class 'dapp'}
\item{tilt.prior}{lofical giving whether the prior should be tilted to mimic an analysis done with a uniform prior on the range(alpha)}
\item{mesh.tilt}{a tuning parameter that controls how exactly tilting is done. Shorter mesh value gives tighter match but will require more Monte Carlo simulations}
\item{nprior}{number of prior draws to be used for display}
\item{ncurves}{number of curves to be shown individually}
\item{simple.layout}{logical indicating if a simpler graphical output should be returned with only predictive visualization}
\item{...}{no addiitonal parameters used at this point}
}
\value{
Gives prior and posterior summaries of the range and average predicted alpha curves
}

\seealso{\code{\link{dapp}}, \code{\link{predict.dapp}} and \code{\link{summary.dapp}}.}

\examples{
\dontrun{
## generate 25 A and 30 B trials with rate functions
##    lambda.A(t) = 160*exp(-2*t/1000) + 40*exp(-0.2*t/1000)
##    lambda.B(t) = 40*exp(-2*t/1000)
## where time t is measured in ms. Then, generate 40 AB trials,
## roughly half with flat weight curves with a constant intensity
## either close to A, or close to B or close to the 50-50 mark,
## (equally likely). The remaining curves are sinusoidal
## that snake between 0.01 and 0.99 with a period randomly
## drawn between 400 and 1000

ntrials <- c(nA=25, nB=30, nAB=40)
flat.range <- list(A=c(0.85, 0.95),
                   B=c(0.05, 0.15),
                   mid=c(0.45,0.55))
flat.mix <- c(A=1/3, B=1/3, mid=1/3)
wavy.span <- c(0.01, 0.99)
wavy.period <- c(400, 1000)

T.horiz <- 1000
rateB <- 40 * exp(-2*(1:T.horiz)/T.horiz)
rateA <- 4*rateB + 40 * exp(-0.2*(1:T.horiz)/T.horiz)

synth.data <- synthesis.dapp(ntrials = ntrials, pr.flat = 0.5,
                             intervals = flat.range, wts = flat.mix,
                             span = wavy.span, period.range = wavy.period,
                             lambda.A=rateA, lambda.B=rateB)

## Visualize data and generated binned spike counts
spike.counts <- mplex.preprocess(synth.data$spiketimes, visualize=TRUE)

## Fit the DAPP model to data
fit.post <- dapp(spike.counts, verbose=FALSE)

## Visualize model fit
plot(fit.post)

## Post process results to assign second order stochasticity labels
summary(fit.post)
}
}
\keyword{programming}
