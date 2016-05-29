/******************************************************************************************
 *
 *  Creates wrapper function to call anclik from R. Likelihood based on 
 *  Alkema, Raftery, Clark Ann Appl Stat 2007. (http://dx.doi.org/10.1214/07-AOAS111)
 *
 *  GPLv3, no warranty, etc...
 *
 *  Created by Jeff Eaton on 2014-11-12 <jeffrey.eaton@imperial.ac.uk>
 *  Updated by Jeff Eaton on 2016-05-29: removed dependency on GSL, use R API instead
 *
 *****************************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>

struct anclik_integrand_param{
  size_t numsites;      // number of sites
  size_t *nobs_site;    // number of observations per site
  double **dst;         // array of vectors: difference between observed and model (transformed) for each site
  double **vst;         // array of vectors: variance for each site (transformed)
  double *x;            // vectors of length maxobs: hold copy of dst for each ldmvnorm calculation
  double *mean;         // vector of 0s, length maxobs
  double *sigma;        // vector of length maxobs^2: memory to create covariance matrix
  double s2_pr_alpha;
  double s2_pr_beta;
};

void anclik_integrand(double *val, int n, void * params);

SEXP anclikR(SEXP s_dst, SEXP s_vst, SEXP s_s2_pr_alpha, SEXP s_s2_pr_beta){
  // s_dst: list of vectors with differences between transformed observed and modelled prevalences
  // s_vst: list of vectors with clinic level variances (transformed) for each site
  // s_s2_pr_alpha: parameter for inverse-gamma prior on variance of clinic-level effects
  // s_s2_pr_beta: parameter for inverse-gamma prior on variance of clinic-level effects

  size_t numsites = length(s_dst);
  size_t *nobs_site = (size_t*) R_alloc(numsites, sizeof(size_t));
  size_t maxobs = 0;

  double **dst = (double**) R_alloc(numsites, sizeof(double*));
  double **vst = (double**) R_alloc(numsites, sizeof(double*));

  for(size_t site = 0; site < numsites; site++){
    nobs_site[site] = length(VECTOR_ELT(s_dst, site));
    dst[site] = REAL(VECTOR_ELT(s_dst, site));
    vst[site] = REAL(VECTOR_ELT(s_vst, site));
    if(maxobs < nobs_site[site])
      maxobs = nobs_site[site];
  }

  double *x = (double*) R_alloc(maxobs, sizeof(double));
  double *mean = (double*) R_alloc(maxobs, sizeof(double));
  double *sigma = (double*) R_alloc(maxobs*maxobs, sizeof(double));
  for(size_t i=0; i < maxobs; i++)
    mean[i] = 0.0;

  struct anclik_integrand_param param = {numsites, nobs_site, dst, vst, x, mean, sigma, *REAL(s_s2_pr_alpha), *REAL(s_s2_pr_beta)};

  // create R object for return and call anclik
  SEXP s_val;
  PROTECT(s_val = allocVector(REALSXP, 1));

  // numerical integration
  // Note: integration limits and parameters taken from Leontine's R implementation.
  double a=1.0e-15, b=0.3, epsabs=0.0001220703, epsrel=0.0001220703, err;
  int neval, ier, last;
  int limit = 1000;
  int lenw = 4 * limit;
  int *iwork = (int *) R_alloc(limit, sizeof(int));
  double *work = (double *) R_alloc(lenw,  sizeof(double));
 
  Rdqags(anclik_integrand, &param, &a, &b, &epsabs, &epsrel,
         REAL(s_val), &err, &neval, &ier, &limit, &lenw, &last, iwork, work);


  UNPROTECT(1);

  return s_val;
}


// Log density of multivariate normal distribution
// NOTE: OVERWRITES *x AND *sigma!!!
double ldmvnorm(int n, double *x, double *mean, double *sigma){
  int info, inc = 1;

  // x = x - mean
  for(int i=0; i<n; i++)
    x[i] = x[i] - mean[i];
  
  F77_NAME(dpotrf)("U", &n, sigma, &n, &info); // chol(sigma);
  F77_NAME(dtrsv)("U", "T", "N", &n, sigma, &n, x, &inc); // inv(L) %*% (x-mean)
  double rss = F77_NAME(ddot)(&n, x, &inc, x, &inc); // (x-mean) %*% inv(sigma) %*% (x-mean)

  double logsqrtdet = 0.0;
  for(int i=0; i<n; i++)
    logsqrtdet += log(sigma[i*n+i]);

  return -logsqrtdet - n*M_LN_SQRT_2PI - 0.5*rss;
}

// likelihood integrand 
void anclik_integrand(double *val, int n, void * params){

  struct anclik_integrand_param *p = (struct anclik_integrand_param *) params;

  for(size_t iter = 0; iter < n; iter++){
    double s2 = val[iter];

    double integrand = 0.0;
    for(size_t i = 0; i < p->numsites; i++){
      size_t nobs = p->nobs_site[i];
      for(size_t j = 0; j < nobs; j++){
	p->x[j] = p->dst[i][j];
	for(size_t k = 0; k < j; k++) // ldmvnorm only need upper-triangle of sigma
	  p->sigma[j*nobs+k] = s2;
	p->sigma[j*nobs+j] = s2+p->vst[i][j];
      }

      integrand += ldmvnorm(nobs, p->x, p->mean, p->sigma);
    }

    integrand += log(s2)*(-p->s2_pr_alpha-1.0) - 1.0/(p->s2_pr_beta * s2);  // Note: omits normalising constant beta^alpha/Gamma(alpha).
    val[iter] = exp(integrand);
  }

  return;
}
