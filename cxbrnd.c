#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include "lehmer.h"


/* Add function to communicate usage information */ 
int usage() {
  printf("CXBrnd is a program to calculate the remaining Cosmic X-ray Background\n");
  printf("flux after the exclusion of point sources up to a certain limit through\n");
  printf("Monte Carlo simulations.\n\n");
  printf("  Usage:\n");
  printf("    ./cxbrnd <flux limit> <area> <niter>\n\n");
  printf("  Example:\n");
  printf("    ./cxbrnd 3.E-15 0.0549 1000\n\n");
  exit(3);
}


/* 
  Program to calculate the variance of background estimates
*/

int main (int argc, char *argv[]) {
  
  double dlim, smin, smax, area;
  double nleh, nsrc, range, s, dns, rej;
  double *x, *y, *tflx, mflx, vflx, low, up, skew;
  unsigned int n;
  int i,j,k,niter;
  
  /* Set GSL random number generator */
  const gsl_rng_type * T;
  gsl_rng * r;
  
  /* Check command line input */
  if (argc != 4) {
    printf("Error: illegal input on command line.");
    usage();
  } 
    
  /* Set basic parameters for the simulations       */
  dlim = atof(argv[1]);  /* Detection limit of telescope   */
  smin = 1.4E-16;        /* Lower limit of log N - log S   */
  smax = 1.0E-13;        /* Maximum flux for log N - log S */  
  area = atof(argv[2]);  /* Field of view in degrees       */     
  
  
  /* Get number of sources from Log N - Log S */
  NsrcLehmer(&smin,&smax,&nleh); 
  nleh = nleh * area;
  
  
  /* Set-up GSL random number generator */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  
  /* Monte Carlo generation of sources */
  niter=atof(argv[3]);
  
  tflx = (double *) malloc(sizeof(double[niter]));
  
  dns = 1.E-20;

  for (i=0;i<niter;i++) {
  
    /* Determine number of sources in area */
    n = gsl_ran_poisson(r,nleh);
    nsrc = (double) n;
    
    /* Set total flux to unresolved fraction */ 
    tflx[i]=3.4e-12*area;
  
    for (j=0;j<n;j++) {
    
      rej = 1.;
    
      while (rej > dns) {
        rej = gsl_rng_uniform(r);
        s = smax * gsl_rng_uniform_pos(r) + smin;
        dns = dNdS(&s);
      };
  
      if (s < dlim) {  
        tflx[i] = tflx[i] + s;
      };
      
    };
  };
  
  /* Close random number generator */
  gsl_rng_free(r);
  
  
  /* Calculate mean and variance */
    
  /* Mean */
  mflx=0.;
  for (i=0;i<niter;i++) {
    mflx = mflx + tflx[i]/((double) niter);
  };
  
  
  /* Variance */
  vflx=0.;
  skew=0.;
  for (i=0;i<niter;i++) {
    vflx = vflx + (pow(tflx[i],2.)-pow(mflx,2.))/((double) niter);
    skew = skew + (pow(tflx[i],3.)-pow(mflx,3.))/((double) niter); 
  };
  
  skew = skew * pow(vflx,-3./2.);
  
  /*
  printf("Read Serr 2\n");
  printf("! Mean, variance and skewness: %e %e %e\n",mflx,sqrt(vflx),skew); 
  */
  printf("%e %e %e\n",dlim,mflx,sqrt(vflx));
  
  
  /* Allocate memory for Histogram */
  gsl_histogram *h = gsl_histogram_alloc(30);
  
  gsl_histogram_set_ranges_uniform(h,mflx-4*sqrt(vflx),mflx+4*sqrt(vflx));
  
  for (i=0;i<niter;i++) {
    gsl_histogram_increment(h,tflx[i]);
  };
  
  x = (double *) malloc(sizeof(double[30]));
  y = (double *) malloc(sizeof(double[30]));
  
  for (j=0;j<30;j++) {
    gsl_histogram_get_range(h,j,&low,&up);
    x[j]=(low+up)/2.;
    y[j]=gsl_histogram_get(h,j);
    /* printf("%e  %f  %f\n",x[j],y[j],sqrt(y[j])); */
  };
  
  free(x);
  free(y);  
  gsl_histogram_free(h);
  
};
