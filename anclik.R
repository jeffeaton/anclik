system("R CMD SHLIB -lgsl -lgslcblas anclikR.c anclik.c mvrandist.c")
dyn.load(paste("anclikR", .Platform$dynlib.ext, sep=""))
library(mvtnorm) # required for calling fnANClik with VERSION = "R"

fnPrepareANCLikelihoodData <- function(anc.prev, anc.n, anchor.year = 1970L, return.data=TRUE){
    ## anc.prev: matrix with one row for each site and column for each year
    ## anc.n: sample size, matrix with one row for each site and column for each year
    ## anchor.year: year in which annual prevalence output start -- to determine index to compare data
    ## NOTE: requires year to be stored in column names of anc.prev

    anc.prev <- anc.prev[apply(!is.na(anc.prev), 1, sum) > 0,] # eliminate records with no observations
    anc.n <- anc.n[apply(!is.na(anc.n), 1, sum) > 0,] # eliminate records with no observations

    ancobs.idx <- mapply(intersect, apply(!is.na(anc.prev), 1, which), apply(!is.na(anc.n), 1, which))  
    ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

    anc.years.lst <- lapply(ancobs.idx, function(i) as.integer(colnames(anc.prev)[i]))
    anc.prev.lst <- setNames(lapply(1:length(ancobs.idx), function(i) as.numeric(anc.prev[i, ancobs.idx[[i]]])), rownames(anc.prev))
    anc.n.lst <- setNames(lapply(1:length(ancobs.idx), function(i) as.numeric(anc.n[i, ancobs.idx[[i]]])), rownames(anc.n))
    
    x.lst <- mapply(function(p, n) (p*n+0.5)/(n+1), anc.prev.lst, anc.n.lst)
    W.lst <- lapply(x.lst, qnorm)
    v.lst <- mapply(function(W, x, n) 2*pi*exp(W^2)*x*(1-x)/n, W.lst, x.lst, anc.n.lst)
    anc.idx.lst <- lapply(anc.years.lst, "-", anchor.year-1)  ## index of observations relative to output prevalence vector


    anclik.dat <- list(W.lst = W.lst,
                       v.lst = v.lst,
                       anc.idx.lst = anc.idx.lst)
    
    if(return.data){ ## Return the data matrices in the list (for convenience)
      anclik.dat$anc.prev <- anc.prev
      anclik.dat$anc.n <- anc.n
    }

    return(anclik.dat)
  }

fnANClik <- function(qM, anclik.dat, s2.pr.alpha = 0.58, s2.pr.beta = 93, VERSION="C"){
    ## qM: vector of probit-transformed annual prevalences (starting in anchor.year specified for anclik.dat)
    ## anclik.dat: list including transformed ANC prevalence data
    ## s2.pr.alpha: parameter for inverse-gamma prior on ANC site-level effects
    ## s2.pr.beta: parameter for inverse-gamma prior on ANC site-level effects
    
    d.lst <- mapply(function(w, idx) w - qM[idx], anclik.dat$W.lst, anclik.dat$anc.idx.lst)

    if(VERSION == "R"){
        V.lst <- lapply(anclik.dat$v.lst, function(x) diag(x, nrow=length(x)))
        return(integrate(Vectorize(function(s2)
                                   exp(sum(mapply(dmvnorm, x=d.lst, sigma = lapply(V.lst, function(m) s2+m), MoreArgs=list(log=TRUE))))*s2^(-s2.pr.alpha-1)*exp(-1/(s2.pr.beta*s2))), 1e-15, 0.3, subdivisions=1000, stop.on.error=FALSE)$value)
    }
    
    return(.Call("anclikR", d.lst, anclik.dat$v.lst, s2.pr.alpha, s2.pr.beta))
}
