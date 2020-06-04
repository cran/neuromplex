\name{synthesis.dapp}
\Rdversion{0.2}
\alias{synthesis.dapp}
\title{Simulate Multiplexing Data for DAPP Analysis}
\description{Simulate spike trains from controlled DAPP setting with flat and sinusoidal weight curves}
\usage{ 
synthesis.dapp(ntrials = c(10, 10, 10), time.bins = 0:1000, lambda.A = 400,
     lambda.B = 100, pr.flat = 0.5, intervals = list(c(0,1)),
     wts = 1, span = c(0,1), period.range = c(400, 1000))
} 
\arguments{
\item{ntrials}{a vector of 3 elements giving the trial counts for conditions A, B and AB}
\item{time.bins}{time bins (in ms) giving the break points of the time bins in which Poisson draws should be made to mimic a Poisson process generation}
\item{lambda.A}{a flat intensity (in Hz) for condition A}
\item{lambda.B}{a flat intensity (in Hz) for condition B}
\item{pr.flat}{proportion of flat weight curves to be generated}
\item{intervals}{a list of sub-intervals (each represented by the 2-vector giving the sub-interval end-points) which determine the ranges of the flat weight curves}
\item{wts}{the relative weights of the sub-intervals above}
\item{span}{a two-vector giving the range of the sinusoidal weight curves}
\item{period.range}{the range from which the sinusoidal periods are drawn randomly (and uniformly)}
}
\value{
Returns a list containting the following items.

\item{spiketimes}{a list with 3 elements giving the 3 sets of spiketimes associated with experimental conditions A, B and AB}
\item{alphas}{true underlying weight curves for each AB trial}
\item{lambdas}{corresponding intensity curves for each AB trial}
\item{time.pts}{time points associated with alphas and lambdas}
}

\examples{
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
}

\keyword{programming}
