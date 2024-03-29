expand.grid <- function(x, y){
  nx <- ifelse(is.null(dim(x)), length(x), nrow(x))
  ny <- ifelse(is.null(dim(y)), length(y), nrow(y))
  cbind(kronecker(rep(1, ny), x), kronecker(y, rep(1, nx)))
}

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))
logdiff <- function(lx) return(ifelse(diff(lx) < 1, lx[1] + log(expm1(diff(lx))), lx[2] + log(1 - exp(lx[1]-lx[2]))))
bin.counter <- function(x, b) return(diff(sapply(b, function(a) sum(x <= a))))

log.pm <- function(x, a, b){
  a.x <- a + sum(x)
  b.x <- b + length(x)
  mu.map <- a.x / b.x
  return(sum(dpois(x, mu.map, log = TRUE)) - diff(dgamma(mu.map, c(a, a.x), c(b, b.x), log = TRUE)))
}

log.pm.ibf <- function(x, a, b) return(mean(sapply(1:length(x), function(j) log.pm(x[-j], x[j] + a, 1 + b))))

log.pm.mix <- function(alloc, z, a, b, c1, c2){
  term1 <- 0; term2 <- 0
  if(!all(alloc == 2)) term1 <- log.pm(z[alloc == 1], a[1], b[1])
  if(!all(alloc == 1)) term2 <- log.pm(z[alloc == 2], a[2], b[2])
  term3 <- lbeta(c1 + sum(alloc == 1), c2 + sum(alloc == 2)) - lbeta(c1, c2)
  return(term1 + term2 + term3)
}

rmult <- function(p) sample(length(p), 1, prob = p)

test.chisq <- function(x, k){
  mu <- mean(x)
  n <- length(x)
  cuts <- qpois((0:k) / k, mu)
  obs.c <- bin.counter(x, cuts)
  exp.c <- n / k
  return(sum((obs.c - exp.c)^2 / exp.c))
}

pval.chisq <- function(x, method=c("variance", "bincount")[1]){
  n <- length(x)
  if(method == 'bincount'){
    k <- max(3, floor(n/5))
    Q.obs <- test.chisq(x,k)
    Q.samp <- replicate(1e4, test.chisq(rpois(length(x), mean(x)),k))
    pval <- mean(Q.samp > Q.obs)
  } else {
    Tstat <- (n-1)*var(x)/mean(x)
    pval <- pchisq(Tstat, df=n-1, lower.tail=FALSE)
  }
  return(pval)
}

poisson.tests <- function(xA, xB, xAB, labels = c("A", "B", "AB"), 
                          remove.zeros = FALSE, gamma.pars = c(0.5, 2e-10), 
                          beta.pars = c(0.5, 0.5), nMC = 1000,
                          plot = FALSE, add.poisson.fits=FALSE, 
                          method.screen = c('variance', 'bincount'),
                          ...){
  
  a <- gamma.pars[1]; b <- gamma.pars[2]
  c1 <- beta.pars[1]; c2 <- beta.pars[2]
  
  if(remove.zeros){
    xA <- xA[xA != 0]
    xB <- xB[xB != 0]
    xAB <- xAB[xAB != 0]
  }
  
  nA <- length(xA)
  nB <- length(xB)
  nAB <- length(xAB)
  if(nA == 0 | nB == 0) stop("not enough data in single sound")
  x.max <- max(max(xA), max(xB), max(xAB))
  
  total.counts <- list(A=xA, B=xB, AB=xAB)
  names(total.counts) <- labels
  x.means <- sapply(total.counts, mean)
  if(plot){
    AB.palette <- c("#E69F0088", "#56B4E9CC", '#1F968BCC')
    tot.df <- bind_rows(lapply(total.counts, function(z) data.frame(Counts=z)), .id = 'Condition')
    tot.df$Condition <- factor(tot.df$Condition, levels=labels)
    g <- ggplot(tot.df, aes(x=Counts, fill=Condition)) + 
      geom_density(col=0) + 
      scale_fill_manual(values=AB.palette) + 
      theme_minimal() + 
      ylab("Density")
    
    if(add.poisson.fits){
      x.range <- layer_scales(g)$x$range$range
      x.grid <- x.range[1]:x.range[2]
      fn.list <- lapply(total.counts, 
                        function(z) data.frame(x=x.grid, pmf=dpois(x.grid,mean(z))))
      
      dfn <- bind_rows(fn.list,.id="Condition")
      dfn$Condition <- factor(dfn$Condition, levels=labels)
      g <- g + 
        geom_line(data=dfn, aes(x=x, y=pmf, col=Condition)) + 
        scale_color_manual(values=AB.palette)
    }
  }
  pvls <- sapply(list(xA, xB, xAB), pval.chisq, method= method.screen[1])
  ## how different are the two pure trials?
  
  two.poi.ibf <- Vectorize(function(i, j) return(log.pm(xA[-i],xA[i]+a,1+b) + log.pm(xB[-j],xB[j]+a,1+b) - log.pm(c(xA[-i],xB[-j]),xA[i]+xB[j]+a,2+b)))
  lbf.pure <- mean(c(outer(1:length(xA), 1:length(xB), two.poi.ibf)))
  
  ## calculate log-marginal prob for the base model
  ## xAB ~ Poi(mu3) with mu3 ~ Ga(a, b)
  ## NB. This is not to be included as a hypothesis.. only for benchmarking
  
  log.marg.ind <- log.pm.ibf(xAB, a, b)
  
  ## calculate log-marginal prob for the mixture model
  ## xAB ~ p * Poi(mu1) + (1 - p) * Poi(mu2)
  ## p ~ Be(c1, c2)
  ## where mu1 and mu2 are the parameters for xA and xB
  
  b.post <- b + c(length(xA), length(xB))
  a.post <- a + c(sum(xA), sum(xB))
  lp1 <- sapply(xAB, log.pm, a = a.post[1], b = b.post[1])
  lp2 <- sapply(xAB, log.pm, a = a.post[2], b = b.post[2])
  lprob <- cbind(lp1, lp2)
  prob <- exp(lprob - apply(lprob, 1, logsum))
  
  
  if(nAB > 10){
    alloc <- replicate(nMC, apply(prob, 1, rmult))
    lpmx <- apply(alloc, 2, log.pm.mix, z = xAB, a = a.post, b = b.post, c1 = c1, c2 = c2)
    lprop <- apply(alloc, 2, function(ix) return(sum(log(diag(prob[,ix])))))
    log.marg.mix.f <- logsum(lpmx - lprop) - log(1e3)
  } else {
    alloc.all <- c(1, 2)
    if(length(xAB) > 1){
      for(i in 1:(length(xAB) - 1)) alloc.all <- expand.grid(alloc.all, c(1, 2))
    }
    alloc.all <- matrix(alloc.all, ncol = length(xAB))
    lpmx.all <- apply(alloc.all, 1, log.pm.mix, z = xAB, a = a.post, b = b.post, c1 = c1, c2 = c2)
    log.marg.mix.f <- logsum(lpmx.all)
  }
  
  log.marg.mix.s.1 <- sapply(xAB, function(z) log.pm.mix(1, z, a.post, b.post, c1, c2))
  log.marg.mix.s.2 <- sapply(xAB, function(z) log.pm.mix(2, z, a.post, b.post, c1, c2))
  log.marg.mix.s <- apply(cbind(log.marg.mix.s.1, log.marg.mix.s.2), 1, logsum)
  log.marg.mix <- log.marg.mix.f - mean(log.marg.mix.s)
  
  ## calculate log-marginal prob for the betweener model
  ## shorties ~ Poi(mu3)
  ## mu3 ~ Ga(a, b) restricted to (mu1, mu2)
  ## where mu1 and mu2 are the parameters for xA and xB
  
  nsamp <- nMC
  mu1.samp <- rgamma(nsamp, a.post[1], b.post[1])
  mu2.samp <- rgamma(nsamp, a.post[2], b.post[2])
  par.samp <- cbind(mu1.samp, mu2.samp)
  log.p <- function(par, z){
    n <- length(z)
    s <- sum(z) 
    rpar <- range(par)
    return(logdiff(pgamma(rpar, s+a, n+b, log.p = TRUE)) - logdiff(pgamma(rpar,a,b,log.p = TRUE)) - (s+a)*log(n+b) + lgamma(s+a) - sum(lgamma(z+1)))
  } 
  log.marg.ave.f <- logsum(apply(par.samp, 1, log.p, z = xAB)) - log(nsamp)
  log.marg.ave.s <- sapply(xAB, function(zz) logsum(apply(par.samp, 1, log.p, z = zz)) - log(nsamp))
  log.marg.ave <- log.marg.ave.f - mean(log.marg.ave.s)
  
  ## calculate log-marginal prob for the "outside" model
  ## shorties ~ Poi(mu3)
  ## mu3 ~ Ga(a, b) restricted to mu3 > max(mu1, mu2) or mu3 < min(mu1, mu2)
  ## where mu1 and mu2 are the parameters for xA and xB
  
  log.pout <- function(par, z){
    n <- length(z)
    s <- sum(z) 
    rpar <- c(max(par), Inf)
    term1 <- logdiff(pgamma(rpar, s+a, n+b, log.p = TRUE)) - logdiff(pgamma(rpar,a,b,log.p = TRUE)) - (s+a)*log(n+b) + lgamma(s+a) - sum(lgamma(z+1))
    rpar <- c(min(par),0)
    term2 <- logdiff(pgamma(rpar, s+a, n+b, log.p = TRUE, lower.tail =  FALSE)) - logdiff(pgamma(rpar,a,b,log.p = TRUE,lower.tail =  FALSE)) - (s+a)*log(n+b) + lgamma(s+a) - sum(lgamma(z+1))
    return(logsum(c(term1, term2)))
  } 
  log.marg.out.f <- logsum(apply(par.samp, 1, log.pout, z = xAB)) - log(nsamp)
  log.marg.out.s <- sapply(xAB, function(zz) logsum(apply(par.samp, 1, log.pout, z = zz)) - log(nsamp))
  log.marg.out <- log.marg.out.f - mean(log.marg.out.s)
  
  ## calculate log-marginal prob for dominant model
  ## shorties ~ Poi(mu1) or Poi(mu2)  
  ## where mu1 and mu2 are the parameters for xA and xB
  log.marg.dom.1 <- log.pm.ibf(xAB, a.post[1], b.post[1]) 
  log.marg.dom.2 <- log.pm.ibf(xAB, a.post[2], b.post[2])
  log.marg.dom <- max(log.marg.dom.1,log.marg.dom.2)
  
  bfs <- exp(c(mixture = log.marg.mix, intermediate = log.marg.ave, outside = log.marg.out, single = log.marg.dom) - log.marg.ind)
  out <- list(separation.logBF = lbf.pure,
              post.prob = bfs/sum(bfs),
              pois.pvalues = pvls[1:2],
              samp.sizes = c(nA, nB, nAB))
  model.probs <- bfs/sum(bfs)
  models <- c("Mixture", "Intermediate", "Outside", "Single")
  win.model <- which.max(model.probs)
  if(plot){
    title.code <- bquote(BayesFactor[A!=B]==.(signif(exp(lbf.pure),1))~"Winning Model"==.(models[win.model])~.(round(100*model.probs[win.model]))~"%")
    g2 <- g + ggtitle(title.code)
    grid.arrange(g2, ...)
  }
  return(out)
}

