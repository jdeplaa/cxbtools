/* Source file for the Lehmer Log N - Log S function */
/* Jelle de Plaa, SRON Utrecht, 18-07-2013           */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lehmer.h"


/*  Parameters for the 2-8 keV Log N - Log S curve */
parameters agn = { 5.7313e+16, 1.32, 2.55, 6.4e-15 };
parameters gal = { 1.10e+14  , 2.29, 0., 1.0 };
parameters star= { 0.64e+14  , 1.79, 0., 1.0 };


/* 
   =========================================================
   Function to calculate the flux emitted by sources between
   slow and shigh, weighted by the Log N - Log S curve.
   Integral | S (dN/dS) dS equals Omega CXB. 
   See Lehmer et al. (2012) Section 3.4. 
   =========================================================
*/ 

int lehmer(double *fmin, double *fmax, double *flux) {  
  
  const double sref = 1.0e-14;  /* Reference flux                          */
  const double flow = 1.4e-16;  /* Detection limit for the 2-8 keV band    */
                                /* See Hickox & Markevitch (2006) Sect. 3.2*/
  
  /* Call error if cut is below the detection limit of Chandra in CDF */
  if ( *fmin < flow ) {
    *fmin = 1.4e-16;
    *flux = 3.4e-12;
  } else {
    *flux = 0.0;
  }  
  
  /* Contribution from unresolved point sources in Chandra DF */
  /* *flux = 3.4e-12;              erg cm^2 s^-1 deg^-2       */
                                /* Hickox & Markevitch (2006) */
                                   
  /* Contribution from AGN */
  if ( *fmin <= agn.fb && *fmax <= agn.fb ) { 
    *flux = *flux + agn.k * pow(sref,agn.b1) * IntLehmer(fmin,fmax,&agn.b1);
  } else if ( *fmin <= agn.fb && *fmax > agn.fb ) {
    *flux = *flux + agn.k * pow(sref,agn.b1) * IntLehmer(fmin,&agn.fb,&agn.b1);
    *flux = *flux + agn.k * pow(agn.fb/sref,agn.b2-agn.b1) * pow(sref,agn.b2) * IntLehmer(&agn.fb,fmax,&agn.b2);
  } else {
    *flux = *flux + agn.k * pow(agn.fb/sref,agn.b2-agn.b1) * pow(sref,agn.b2) * IntLehmer(fmin,fmax,&agn.b2);
  };

  /* Contribution from Galaxies */
  *flux = *flux + gal.k * pow(sref,gal.b1) * IntLehmer(fmin,fmax,&gal.b1);
  
  /* Contribution from Stars */
  *flux = *flux + star.k * pow(sref,star.b1) * IntLehmer(fmin,fmax,&star.b1);
  
  return 0;
};

/* Integration result of the integral | S^(1-b) dS over slow to shigh */

double IntLehmer(const double *slow, double *shigh, double *beta) 
{
  double result;
  
  result = pow(2.0-*beta,-1.0) * pow(*shigh,2.0-*beta);
  result = result - pow(2.0-*beta,-1.0) * pow(*slow,2.0-*beta);
  
  return result;
};



/* 
   ====================================================
   Function to calculate the number of sources between
   slow and shigh, weighted by the Log N - Log S curve.
   Integral | (dN/dS) dS. See Lehmer et al. (2012) 
   ====================================================
*/ 

int NsrcLehmer(double *fmin, double *fmax, double *nsrc)
{
  
  const double sref = 1.0e-14;  /* Reference flux                          */
  const double flow = 1.4e-16;  /* Detection limit for the 2-8 keV band    */
                                /* See Hickox & Markevitch (2006) Sect. 3.2*/
  
  /* Call error if cut is below the detection limit of Chandra in CDF */
  if ( *fmin < flow ) {
    printf("Error: Flux limit is below Chandra detection limit of 1.4e-16\n");
    exit(1);
  } else {
    *nsrc = 0.0;
  };  
  
  /* Contribution from AGN */
  if ( *fmin <= agn.fb && *fmax <= agn.fb ) { 
    *nsrc = *nsrc + agn.k * pow(sref,agn.b1) * IntNsrc(fmin,fmax,&agn.b1);
  } else if ( *fmin <= agn.fb && *fmax > agn.fb ) {
    *nsrc = *nsrc + agn.k * pow(sref,agn.b1) * IntNsrc(fmin,&agn.fb,&agn.b1);
    *nsrc = *nsrc + agn.k * pow(agn.fb/sref,agn.b2-agn.b1) * pow(sref,agn.b2) * IntNsrc(&agn.fb,fmax,&agn.b2);
  } else {
    *nsrc = *nsrc + agn.k * pow(agn.fb/sref,agn.b2-agn.b1) * pow(sref,agn.b2) * IntNsrc(fmin,fmax,&agn.b2);
  };

  /* Contribution from Galaxies */
  *nsrc = *nsrc + gal.k * pow(sref,gal.b1) * IntNsrc(fmin,fmax,&gal.b1);
  
  /* Contribution from Stars */
  *nsrc = *nsrc + star.k * pow(sref,star.b1) * IntNsrc(fmin,fmax,&star.b1);
  
  return 0;
};

double IntNsrc(const double *slow, double *shigh, double *beta)
{
  double result;
  
  result = pow(1.0-*beta,-1.0) * pow(*shigh,1.0-*beta);
  result = result - pow(1.0-*beta,-1.0) * pow(*slow,1.0-*beta);
  
  return result;
};


/* 
   ====================================================
   Function that returns the dN/dS to 
   find a source in flux interval dS
   ====================================================
*/ 

double dNdS(double *s) 
{
  double dns;
  const double sref = 1.0e-14;  /* Reference flux */
  
  if ( *s <= agn.fb ) {
    dns = agn.k * pow(*s/sref,-agn.b1);
  } else { 
    dns = agn.k * pow(agn.fb/sref,agn.b2-agn.b1) * pow(*s/sref,-agn.b2);
  };
  
  dns = dns + star.k * pow(*s/sref,-star.b1);
  
  dns = dns + gal.k * pow(*s/sref,-gal.b1);
  
  dns = dns / 1.E+20 ;
  
  return dns;
};

