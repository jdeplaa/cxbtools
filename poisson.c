#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Poisson function using rejection sampling
   This function draws a random Poisson deviate with
   expectation value mu. For expectation values 
   greater than 1.E+7 we use the normal distribution */
      
int poisson(double mu) {
  const gsl_rng_type * T;
  gsl_rng * r;
  
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  /* unsigned int gsl_ran_poisson (const gsl_rng * r, double mu) */
  
  unsigned int k = gsl_ran_poisson(r,mu);
  
  gsl_rng_free (r);
  return k;
  
}
