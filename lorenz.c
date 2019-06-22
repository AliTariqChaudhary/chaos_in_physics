#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Author: Ali Chaudhary ID: 67694450
  Physics 53

  Description: This code seeks to model a Lorenz attractor.  
               The Lorenz attractor is a Chaotic solution to 
               the Lorenz system. We will use it to demonstrate
               chaos.
   
 Relevent Equations: 

 dx/dt = sigma(y-x)

 dy/dt = x(rho -z) - y

 dz/dt = xy-Beta*z

*/

//System of 3 differential equations
//dx/dt
double fx(double sigma, double x, double y)
{
  return sigma*(y-x);
}
//dy/dt
double fy(double rho, double x, double y, double z)
{
  return x*(rho-z)-y;
}
double fz(double beta, double x, double y, double z)
{
  return x*y-beta*z;
}

void rk4(double ti, double xi, double yi, double zi,
	 double sigma, double rho, double beta,
	 double tf, double* xf, double* yf, double* zf)
{
  double k1x, k2x, k3x, k4x;
  double k1y, k2y, k3y, k4y;
  double k1z, k2z, k3z, k4z;

  double h = tf - ti;
  double t = ti;
  t+=0;
  
  //Do it for each 3 equations
  k1x = fx(sigma, xi, yi);
  k1y = fy(rho, xi, yi, zi);
  k1z = fz(beta, xi, yi, zi);

  k2x = fx(sigma, xi+k1x*h/2.0, yi+k1y*h/2.0);
  k2y = fy(rho, xi+k1x*h/2.0, yi+k1y*h/2.0, zi+k1z*h/2.0);
  k2z = fz(beta,xi+k1x*h/2.0, yi+k1y*h/2.0, zi+k1z*h/2.0);

  k3x = fx(sigma, xi+k2x*h/2.0, yi+k2y*h/2.0);
  k3y = fy(rho, xi+k2x*h/2.0, yi+k2y*h/2.0, zi+k2z*h/2.0);
  k3z = fz(beta,xi+k2x*h/2.0, yi+k2y*h/2.0, zi+k2z*h/2.0);

  k4x = fx(sigma, xi+k3x*h, yi+k3y*h/2.0);
  k4y = fy(rho, xi+k3x*h, yi+k3y*h, zi+k3z*h);
  k4z = fz(beta,xi+k3x*h, yi+k3y*h, zi+k3z*h);

  *xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)* h/6.0;
  *yf = yi + (k1y + 2.0*(k2y+k3y) + k4y)* h/6.0;
  *zf = zi + (k1z + 2.0*(k2z+k3z) + k4z)* h/6.0;

}



int main()
{
  //Firstly I intend to explain to the user what this program does.
  printf("\n\n ----------------------------------------------\n");
  printf("|                                              |\n");
  printf("|      WE WILL USING A LORENZ ATTRACTOR TO     |\n");
  printf("|      DEMONSTRATE CHAOS IN PHYSICS            |\n");
  printf("|                                              |\n");
  printf(" ----------------------------------------------\n\n\n"); 
  
  double sig = -1;
  while(sig<0)
    {
      printf("\nPlease enter a value for sigma greater than zero:");
      scanf("%lf",&sig);
      if(sig <0)
	printf("That is less than zero, try again");
    }
  double rho = -1;
  while(rho<0)
    {
      printf("\nPlease enter a value for rho greater than zero:");
      scanf("%lf",&rho);
      if(rho <0)
	printf("That is less than zero, try again");
    }
  double beta = -1;
  while(beta<0)
    {
      printf("\nPlease enter a value for beta greater than zero:");
      scanf("%lf",&beta);
      if(beta <0)
	printf("That is less than zero, try again");
    }
  double Nmax = -1;
  while(Nmax<0)
    {
      printf("\nHow many seconds do you want to plot? (100 is usually a good number:");
      scanf("%lf",&Nmax);
      if(Nmax <0)
	printf("That is less than zero, try again");
    }




  //Now for the implementation of the RK4 method
   double xi,yi,zi,xf,yf,zf;
   double dt = 0.01;
   double tf;



   xi = 0.1;
   yi = 0;
   zi = 0;
   double ti = 0;
   tf = ti;
   xf = xi;
   yf =yi;
   zf = zi;

   //Set up to print to file
   FILE* fPtr;
   char fileName[] = "lorenz.dat";
   fPtr = fopen(fileName,"w");

   while(ti<Nmax)
     {
       tf+= dt;
       rk4(ti, xi, yi, zi, sig, rho, beta, tf, &xf, &yf, &zf);
       fprintf(fPtr, "%lf  %lf  %lf\n", xf, yf, zf);
       ti = tf;
       xi = xf;
       yi = yf;
       zi = zf;
     }

   fclose(fPtr);
   printf("THE DATA HAS BEEN SAVED TO 'lorenz.dat'!\n");
   return 0;
}
