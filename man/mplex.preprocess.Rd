\name{mplex.preprocess}
\Rdversion{0.2}
\alias{mplex.preprocess}
\title{Preprocessing Neural Multiplexing Data}
\description{Preprocess nueral spike train recording to preapre binned spike counts suitable for DAPP analysis}
\usage{ 
mplex.preprocess(spiketimes, start.time=0, end.time=1e3, bw=50,
    remove.zeros=FALSE, visualize=TRUE, ...)
} 
\arguments{
\item{spiketimes}{a list with 3 elements giving the 3 sets of spiketimes associated with experimental conditions A, B and AB}
\item{start.time}{starting time for the observation window. See details below}
\item{end.time}{ending time of the observations window. See details below}
\item{bw}{bin width (in ms) used for binning. A single bin is used when bw equals or exceeds the length of the observation period (end.time - start.time). Single bin analysis is same as total spike count analysis}
\item{remove.zeros}{logical indicating if trials with zero spike counts should be removed from the analysis}
\item{visualize}{logical indicating if a graphical summary should be produced to visualize the three sets of trials}
\item{...}{additional commands to be passed on to grid.arrange() for plotting. For example, adding 'top="PLOT TITLE"' will add a title at the top of the combined plot. See \code{grid.arrange} for more details.}
}

\value{
Returns a list containting the following items.

\item{Acounts}{binned spike counts under condition A presented as a matrix. Rows are bins, columns are replicates (trials). In case of single bin analysis, i.e., with bw equal or larger than total observation window length, a vector of counts is returned.}
\item{Bcount}{binned spike counts under condition B}
\item{ABcounts}{binned spike counts under condition AB}
\item{bin.mids}{an array giving the mid-points of the time bins}
\item{bin.width}{a scalar giving the bin width}
\item{time.horizon}{a vector of length 2 giving the start and the end times of the observation period}
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

## Visualize data and generate binned spike counts
spike.counts <- mplex.preprocess(synth.data$spiketimes, visualize=TRUE, 
top="Synthetic data: bin size=50ms")

\dontrun{
## Visualize total spike counts data
spike.counts <- mplex.preprocess(synth.data$spiketimes, bw=Inf, visualize=TRUE, 
top="Synthetic data: total spike counts")
}
}

\keyword{programming}
