fnPrepareANCLikelihoodData <- function(anc.prev, anc.n, anc.used = TRUE, anchor.year = 1970L, return.data=TRUE){
    ## anc.prev: matrix with one row for each site and column for each year
    ## anc.n: sample size, matrix with one row for each site and column for each year
    ## anchor.year: year in which annual prevalence output start -- to determine index to compare data
    ## NOTE: requires year to be stored in column names of anc.prev

    anc.prev <- anc.prev[anc.used,,drop=FALSE]  # keep only used sites
    anc.n <- anc.n[anc.used,,drop=FALSE]        # keep only used sites
  
    anc.prev <- anc.prev[apply(!is.na(anc.prev), 1, sum) > 0,,drop=FALSE] # eliminate records with no observations
    anc.n <- anc.n[apply(!is.na(anc.n), 1, sum) > 0,,drop=FALSE] # eliminate records with no observations

    ancobs.idx <- mapply(intersect, lapply(as.data.frame(t(!is.na(anc.prev))), which),
                         lapply(as.data.frame(t(!is.na(anc.n))), which), SIMPLIFY=FALSE)
    ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

    anc.years.lst <- lapply(ancobs.idx, function(i) as.integer(colnames(anc.prev)[i]))
    anc.prev.lst <- setNames(lapply(1:length(ancobs.idx), function(i) as.numeric(anc.prev[i, ancobs.idx[[i]]])), rownames(anc.prev))
    anc.n.lst <- setNames(lapply(1:length(ancobs.idx), function(i) as.numeric(anc.n[i, ancobs.idx[[i]]])), rownames(anc.n))
    
    x.lst <- mapply(function(p, n) (p*n+0.5)/(n+1), anc.prev.lst, anc.n.lst, SIMPLIFY=FALSE)
    W.lst <- lapply(x.lst, qnorm)
    v.lst <- mapply(function(W, x, n) 2*pi*exp(W^2)*x*(1-x)/n, W.lst, x.lst, anc.n.lst, SIMPLIFY=FALSE)
    anc.idx.lst <- lapply(anc.years.lst, "-", anchor.year-1)  ## index of observations relative to output prevalence vector


    anclik.dat <- list(W.lst = W.lst,
                       v.lst = v.lst,
                       n.lst = anc.n.lst,
                       anc.idx.lst = anc.idx.lst)
    
    if(return.data){ ## Return the data matrices in the list (for convenience)
      anclik.dat$anc.prev <- anc.prev
      anclik.dat$anc.n <- anc.n
    }

    return(anclik.dat)
  }

fnANClik <- function(qM, anclik.dat, v.infl=0, s2.pr.alpha = 0.58, s2.pr.beta = 93, VERSION="C"){
    ## qM: vector of probit-transformed annual prevalences (starting in anchor.year specified for anclik.dat)
    ## anclik.dat: list including transformed ANC prevalence data
    ## v.infl: additional variance for ANC prevalence observation
    ## s2.pr.alpha: parameter for inverse-gamma prior on ANC site-level effects
    ## s2.pr.beta: parameter for inverse-gamma prior on ANC site-level effects
    
    d.lst <- mapply(function(w, idx) w - qM[idx], anclik.dat$W.lst, anclik.dat$anc.idx.lst, SIMPLIFY=FALSE)
    v.lst <- lapply(anclik.dat$v.lst, "+", v.infl)

  anc_resid_lik(d.lst, v.lst, s2.pr.alpha, s2.pr.beta, VERSION)
}

#' Likelihood for integrated ANC likelihood on residuals. (Alkema, Raftery, Clark equation 10)
#'
#' @param d.lst residuals of observed ANC prevalence compared to predicted
#' @param v.lst approximated variance for each ANC observation
#' @param s2.pr.alpha parameter for inverse-gamma prior on ANC site-level effects
#' @param s2.pr.beta parameter for inverse-gamma prior on ANC site-level effects
#' @param VERSION flag for evaluating C or R implementation of likelihood (for debugging)
#' @return integrated ANC likelihood (double)
anc_resid_lik <- function(d.lst, v.lst, s2.pr.alpha=0.58, s2.pr.beta=93, VERSION="C"){

  if(VERSION == "R"){
     if (!requireNamespace("mvtnorm", quietly = TRUE))
       stop("Package mvtnorm needed to call R version of fnANClik.", call. = FALSE)
     V.lst <- lapply(v.lst, function(x) diag(x, nrow=length(x)))
     return(integrate(Vectorize(function(s2)
       exp(sum(mapply(mvtnorm::dmvnorm, x=d.lst, sigma = lapply(V.lst, function(m) s2+m), MoreArgs=list(log=TRUE))))*s2^(-s2.pr.alpha-1)*exp(-1/(s2.pr.beta*s2))), 1e-15, 0.3, subdivisions=1000, stop.on.error=FALSE)$value)
  }
  
  return(.Call("anclikR", d.lst, v.lst, s2.pr.alpha, s2.pr.beta, PACKAGE="anclik"))
}
  



sample.b.one <- function(d, v, s2.pr.alpha = 0.58, s2.pr.beta = 93){
  ## Use rejection sampling to sample clinic level effect (Alkema, Raftery, Clark 2007)
  ## p(b.s | M, W.s) \propto N( d.st, v.st) * (b.s^2/2 + 1/beta2)^(-alpha-1/2)
  ## 1) sample from normal distribution with weighted mean and variance
  ## 2) reject based on second term of product
  
  max.val <- (1/s2.pr.beta)^(-s2.pr.alpha-0.5)  # maximized when b=0, to normalize rejection sampling
  b <- Inf
  while(runif(1) > (0.5*b^2 + 1/s2.pr.beta)^(-s2.pr.alpha-0.5) / max.val)
    b <- rnorm(1, sum(d/v)/sum(1/v), 1/sum(1/v))
  return(b)
}

sample.b.site <- function(qM, anclik.dat, s2.pr.alpha = 0.58, s2.pr.beta = 93){
  ## Sample b.s values for all clinics
  ## parameters defined the same as fnANClik
  
  d.lst <- mapply(function(w, idx) w - qM[idx], anclik.dat$W.lst, anclik.dat$anc.idx.lst, SIMPLIFY=FALSE)
  return(mapply(sample.b.one, d.lst, anclik.dat$v.lst, s2.pr.alpha, s2.pr.beta))
}

sample.sigma2 <- function(b.site, s2.pr.alpha = 0.58, s2.pr.beta = 93){
  ## Sample from posterior for variance of clinic random effects, conditional
  ## on sample of random effect values.
  
  ## b.site = matrix of b.site values (S x nsample)
  ## Note: scale parameterization, b0 = 1/s2.pr.beta
  ## p(sigma2 | b[1:S], a0, b0) ~ InvGamma(a0 + S/2, b0 + sum(b^2)/2)

  if(is.vector(b.site))
    b.site <- matrix(b.site)
  nsample <- ncol(b.site)
  nsites <- nrow(b.site)
  sigma2 <- 1/rgamma(nsample, shape = s2.pr.alpha+nsites/2, rate = 1/s2.pr.beta+colSums(b.site^2)/2)
  return(sigma2)
}


sample.pred.site <- function(qM, b.site, anclik.dat, v.infl=0){
  ## Sample predicted prevalences in the same years as observed ANC prevalences in each site

  ## qM: vector of probit-transformed annual prevalences (starting in anchor.year specified for anclik.dat)
  ## b.site: vector of site-level random effects
  ## v.infl: additive variance term to (tranformed) binomial variance
  ## Note: b.site must be drawn from same posterior sample as qM (b.site | qM)!!!

  ## site-level fitted values
  fit.site <- mapply(function(b, idx) qM[idx] + b, b.site, anclik.dat$anc.idx.lst, SIMPLIFY=FALSE)

  ## variance fitted values
  vpred.site <- mapply(function(pred, n) 2*pi*exp(pred^2)*pnorm(pred)*(1-pnorm(pred)) / n + v.infl, fit.site, anclik.dat$n.lst, SIMPLIFY=FALSE)

  ## predicted values (probit scale)
  pred.site <- mapply(rnorm, sapply(fit.site, length), fit.site, lapply(vpred.site, sqrt), SIMPLIFY=FALSE)
  
  return(lapply(pred.site, pnorm)) # natural scale
}
