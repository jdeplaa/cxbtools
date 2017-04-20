/* Include file for the Lehmer Log N - Log S function */
/* Jelle de Plaa, SRON Utrecht, 18-07-2013            */

#ifndef lehmer_h
#define lehmer_h

/*
Define structures to contain the parameters for the Lehmer et al. (2012)
fit to the Log N - Log S curve. See Table 1 of Lehmer et al. (2012).
*/

struct paramdef;

typedef struct paramdef {
  /*
  k      K Normalisation for the power law
  b1     Beta 1 power law index below the break
  b2     Beta 2 power law index above the break
  fb     Flux where the break occurs.
  */
  double k, b1, b2, fb;
} parameters;


/*  Parameters for the 2-8 keV Log N - Log S curve */
extern parameters agn ;
extern parameters gal ;
extern parameters star; 


/* Declaration of the function that calculates the integral 
   in section 3.4 of Lehmer et al. (2012) */
    
int lehmer (double *fmin, double *fmax, double *flux);

double IntLehmer(const double *slow, double *shigh, double *beta);

int NsrcLehmer(double *fmin, double *fmax, double *nsrc);

double IntNsrc(const double *slow, double *shigh, double *beta);
 
double dNdS(double *s);  

#endif    
