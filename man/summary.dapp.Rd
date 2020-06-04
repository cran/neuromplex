\name{summary.dapp}
\Rdversion{0.5}
\alias{summary.dapp}
\title{Summary Method for Dynamic Admixture of Poisson Process}
\description{Presents post-processing labels from a DAPP model fit to binned spiking data}
\usage{ 
\method{summary}{dapp}(object, flat.cut = 0.15, wavy.cut = 0.85,
    extreme.cut = 0.25, ...)
} 
\arguments{
\item{object}{a fitted model of the class 'dapp'}
\item{flat.cut}{maximum range allowed to be labelled 'flat'}
\item{wavy.cut}{minimum range allowed to be labelled 'wavy'}
\item{extreme.cut}{for flat curves, maximum deviation from extremes (0 or 1) allowed to be labelled flat.B or flat.A (respectivel)}
\item{...}{additional parameters passed on to the call of \code{\link{predict.dapp}}}
}
\value{
Gives prior and posterior summaries of the range and average predicted alpha curves
}

\details{The summary function analyzes the prior and posterior predictive draws of the weight curves alpha(t). Each draw is assigned with one of the following labels: 'flat.A', 'flat.B', 'flat.Mid', 'wavy', or 'others'. The proportions of these categories are printed for the prior and posterior sets. Additionally, posterior
    draws of alpha(t), for each recorded AB trial, are also analyzed in the same way
    to produce similar labels for each trial, and, the trial is given the label
    that has the maximum posterior probability.
}

\seealso{\code{\link{dapp}}, \code{\link{plot.dapp}} and \code{\link{predict.dapp}}.}

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
