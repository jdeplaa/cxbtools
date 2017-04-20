#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpgplot.h>
#include "lehmer.h"


double eef(double r) {
  double f, noem, tel;
  const double rc = 4.18;
  const double alfa = 1.4;

  tel = 1. - 1./(pow(1.+pow(r/rc,2.),alfa-1.0));
  noem= 1. - 1./(pow(1.+pow(300./rc,2.),alfa-1.0));
  
  f=tel/noem;

  return f;
};

/* Add function to communicate usage information */ 
int usage() {
  printf("CXBopt is a program to calculate optimal source extraction radius\n");
  printf("and flux cut for a certain annulus.\n\n"); 
  printf("  Usage:\n");
  printf("    ./cxbopt <source cnts> <backg counts> <region>\n\n");
  printf("  Example:\n");
  printf("    ./cxbopt 10000. 15000. 0.0549\n\n");
  printf("The <region> area is in square degrees.\n\n");
  exit(3);
}



/* Main Program */

int main (int argc, char *argv[] ) {  
  
  /* Parameters possibly modified by users */ 
  /* double atot = 6.11e-3;     /* Area of the extraction region degree^2 */
  /* double csrc = 70000.;      /* source counts */
  /* double ibkg = 10913.;      /* instrumental background */
  
  double atot, csrc, ibkg;
  double rs;            /* point source extraction radius (arcsec) */
  double slim;          /* Flux limit for excision of point sources (erg/cm2/s) */
  double nsrc, tmp, flux1, flux2, rmax, fmax;
  const double pi = 3.14159265358979323846;
  double snull = 1.4E-16;
  double smax  = 1.E-11;
  int i,j,k, nc;
  float snr[15000]; 
  float c[20], max, min;
  const float tr[6] = {0.5,0.5,0.0,-15.5,0.0,0.03};
  const int   dim[2] = {149,99};
  
  /* PGIMAG color scale */
  const float rl[9] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
  const float rr[9] = {0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  const float rg[9] = {0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  const float rb[9] = {0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
  
  
  if (argc == 4) {
    csrc = atof(argv[1]);
    ibkg = atof(argv[2]);
    atot = atof(argv[3]);
  } else {
    usage();
  }
  
  min=1e+20;
  max=0.0;
  
  for (j=0;j<=dim[0];j++)
  {
    slim = pow(10,tr[3]+tr[5]*j);
    
    NsrcLehmer(&slim,&smax,&nsrc);
    lehmer(&snull,&slim,&flux1);
    lehmer(&slim,&smax,&flux2);
  
    flux1 = flux1 * atot * 1.E+5 * 1.e+3 / 6.408E-9 ;
    flux2 = flux2 * atot * 1.E+5 * 1.e+3 / 6.408E-9 ;
  
    for (i=0;i<=dim[1];i++)
    {
      rs = (tr[0] + tr[1]*i);
      tmp = sqrt(atot*(1. - nsrc*pi*pow(rs/3600.,2.0))) * csrc;
      tmp = tmp / sqrt(csrc + ibkg + flux1 + (1.-eef(rs))*flux2 );  
      k = i + j*100;
      snr[k] = tmp;
      if (snr[k] > max) {
        max=snr[k];
	rmax=rs;
	fmax=slim;
      };
	
      if (snr[k] < min) min=snr[k];
    };
  };
  
  
  /* Typical Nsrc value for this annulus */
  snull=1.E-15;
  slim=3.E-15;
  
  NsrcLehmer(&snull,&slim,&nsrc);
  
  printf("Amount of sources in FOV: %f\n", nsrc*atot); 
  printf("Optimal radius:   %f\n", rmax);
  printf("Optimal flux cut: %e\n",fmax);
  
  /* Plot contours using PGPLOT */
  
  /* Contour levels */
  for (i=0; i<20; i++) 
  {
    c[i]=min + 0.05*(max-min)*i;
  } 
  
  if (cpgbeg(0,"?",1,1) != 1) return EXIT_FAILURE;
  
  /* Set line width to 3 */
  cpgslw(3.0);
  
  /* Draw the box aorund the plot */
  cpgsvp(0.1, 0.9, 0.1, 0.9);
  cpgswin(1.0, 50., tr[3], tr[3]+tr[5]*dim[0] );
  cpgbox("bcnts", 0.0, 0, "bcntsv", 1.0, 10);

  /* Draw image */
  cpgctab(rl,rr,rg,rb,9,1.0,0.5);
  cpgimag(snr,100,150,1,100,1,150,min,max,tr);
  cpgbox("bcnts", 0.0, 0, "bcntsv", 1.0, 10);

  /* Draw contours with labels */
  cpgcont(snr,100,150,1,100,1,150,c,20,tr);  
  
  /*
  cpgconl(snr,100,150,1,100,1,150,c[2],tr,"2",30,10);
  cpgconl(snr,100,150,1,100,1,150,c[5],tr,"5",30,10);
  cpgconl(snr,100,150,1,100,1,150,c[8],tr,"8",30,10);
  */
  
  cpgwedg("RI",2.0,4.0,min,max,"SNR");
  
  /* Draw labels at the axes */  
  cpgsch(1.2);
  cpgmtxt("B",2.8,0.5,0.5,"Point source extraction radius (arcsec)");
  cpgmtxt("L",2.8,0.5,0.5,"Flux cut (erg/cm2/s)");  
  cpgmtxt("T",1.5,0.5,0.5,""); 
   
    
  cpgend();
  
  return 0;
};


