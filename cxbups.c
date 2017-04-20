#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lehmer.h"

/* Add function to communicate usage information */ 
int usage() {
  printf("CXBups is a program to calculate the remaining Cosmic X-ray Background\n");
  printf("flux after the exclusion of point sources up to a certain limit.\n\n"); 
  printf("  Usage:\n");
  printf("    ./cxbups <flux limit> <area>\n\n");
  printf("  Example:\n");
  printf("    ./cxbups 3.E-15 0.0549\n\n");
  exit(3);
}


/* Main Program */

int main (int argc, char *argv[] ) {  
  double flim, fmin, nsrc, var, area, rflux;
  double flux = 0.0;
  
  if (argc == 3) {
    flim = atof(argv[1]);
    area = atof(argv[2]);
    if (flim == 0.0) usage();
  } else {
    usage();
  } 
  
  fmin=0.;
  lehmer(&fmin,&flim,&flux);
  flux = flux * area;
  
  /* Estimate variance on flux assuming Poisson distribution */
  fmin=1.4E-16;
  lehmer(&fmin,&flim,&rflux);
  rflux = rflux * area;
  NsrcLehmer(&fmin,&flim,&nsrc);
  nsrc = nsrc * area;
  var=sqrt(nsrc)*rflux/(nsrc);
  
  /* Convert to SI units for SPEX convenience */
  
  /* printf("Using an input limit of:  %e         erg/cm2/s\n", flim); */
  printf("2-8 keV Flux (cgs):  %e +/- %e  erg/cm2/s \n", flux, var);
  flux = flux * 1.0E-3;
  var  = var * 1.0E-3;
  printf("2-8 keV Flux (SI) :  %e +/- %e  W/m2 \n", flux, var);

  return 0;
};

