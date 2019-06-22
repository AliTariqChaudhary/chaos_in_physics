#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI acos(-1.0)
#define GRAV 980.665 //in cm/s^2 
/* Author: Ali Chaudhary ID: 67694450
  Physics 53

  Description: This code seeks to model a simple pendulum and
               compare its predictable Simple Harmonic Motion
               with the chaotic motion of a double pendulum.
         
               *Since we want to demonstrate the chaos of the
                double pendulum compared to the single pendulum,
                no dampening/air resistence will be used!

               **The equations of motion for a double pendulum
                 actually can't be solved analytically, so the 
                 numerical approximation we will do is the only 
                 way to actually solve it interestingly enough. 
  
  Relevent Equations: 

  Simple Pendulum Equation: d^2(theta)/dt^2 + (g/l)*sin(theta) = 0 
  We will cast this equation to use the Runge Kutta method as:
  d(theta)/dt = v     where v is the angularvelocity
  dv/dt = -g/l *sin(theta)

  For the Double Pendulum we will use similar casting:
 d(theta)/dt = angular velocity

*/


//This is just a function I need to convert an angle and a length to an X cordinate for gnuplot
double X(double angle, double len)
{
  return sin(angle)*len;
}
//Now smame for y
double Y(double angle, double len)
{
  return -cos(angle)*len;
}

//The first order derivative of theta with respect to time for the single pendulum
//It is simply the angular velocity
double single1( double angularvelocity)
{
  return angularvelocity;
}

//Now the derivative of velocity with respect to time
double single2( double len,double theta)
{
  return (-1.0)*GRAV/len*sin(theta);
}


//Now we need the equations for the double pendulum
//There will be 4, 2 for each pendulum

//First pendulum first derivative
double double1a(double angularvelocity1)
{
  return angularvelocity1;
}
//Second pendulum second derivative
double double1v(double t1, double t2, double av1, double av2,
	       double l1, double l2, double m1, double m2)
{
  //t1 is angle one t2 is angle two
  //av is angular velocity
  //l is length
  //m is mass

  double d = t2-t1;//doing this makes calculating the eqn easier
  double answer = m2*l1*av1*av1*sin(d)*cos(d)+m2*GRAV*sin(t2)*cos(d);
  answer += m2*l2*av2*av2*sin(d)-(m1+m2)*GRAV*sin(t1);
  answer /= ((m1+m2)*l1-m2*l1*cos(d)*cos(d));
  return answer;
}

//Second pendulum first derivative
double double2a(double angularvelocity2)
{
  return angularvelocity2;
}

//Second pendulum second derivative
double double2v(double t1, double t2, double av1, double av2,
	       double l1, double l2, double m1, double m2)
{
  //t1 is angle one t2 is angle two
  //av is angular velocity
  //l is length
  //m is mass

  double d = t2-t1;//doing this makes calculating the eqn easier
  double answer = -m2*l2*av2*av2*sin(d)*cos(d)+(m1+m2)*(GRAV*sin(t1)*cos(d)-l1*av1*av1*sin(d)-GRAV*sin(t2));
  answer /= (m1+m2)*l2-m2*l2*cos(d)*cos(d);
  return answer;

}


//The implementation of the Runge Kutta 4 method for the single pendulum
void rk4(double ti, double ai, double vi,
	 double tf, double* af, double* vf, double len)
{
  double tlen = len;
  double k1a, k2a, k3a, k4a;
  double k1v, k2v, k3v, k4v;

  double h = tf - ti;
  double t = ti;
  t+= 0;

  k1a = single1( vi);
  k1v = single2( tlen, ai);

  k2a = single1( vi+k1v*h/2.0);
  k2v = single2( tlen, ai+k1a*h/2.0);

  k3a = single1( vi+k2v*h/2.0);
  k3v = single2( tlen, ai+k2a*h/2.0);

  k4a = single1(vi+k3v*h);
  k4v = single2(tlen, ai+k3a*h);

 
  *af = ai + (k1a + 2.0*(k2a+k3a) + k4a)* h/6.0;
  *vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)* h/6.0;

}

//The implementation of the Rung Kutta 4 method for the double pendulum
void doublerk4(double ti, double ai1, double vi1, double ai2, double vi2, double tf, 
	       double* af1, double* vf1, double* af2, double* vf2, 
	        double l1, double l2, double m1, double m2)
{
  double k1a, k2a, k3a, k4a;
  double k1v, k2v, k3v, k4v;

  double h = tf - ti;
  double t = ti;
  t+= 0;

  //Do one pendulum first then the other
  k1a = double1a( vi1);
  k1v = double1v( ai1, ai2, vi1, vi2, l1, l2, m1, m2);

  k2a = double1a( vi1+k1v*h/2.0);
  k2v = double1v( ai1+k1a*h/2.0, ai2, vi1+k1v*h/2.0, vi2, l1, l2, m1, m2);

  k3a = double1a( vi1+k2v*h/2.0);
  k3v = double1v( ai1+k2a*h/2.0, ai2, vi1+k2v*h/2.0, vi2, l1, l2, m1, m2 );

  k4a = double1a(vi1+k3v*h);
  k4v = double1v(ai1+k3a*h, ai2, vi1+k3v*h, vi2, l1, l2, m1, m2);

 
  *af1 = ai1 + (k1a + 2.0*(k2a+k3a) + k4a)* h/6.0;
  *vf1 = vi1 + (k1v + 2.0*(k2v+k3v) + k4v)* h/6.0;

  //Now for the second one
  k1a = double2a( vi2);
  k1v = double2v( ai1, ai2, vi1, vi2, l1, l2, m1, m2);

  k2a = double2a( vi2+k1v*h/2.0);
  k2v = double2v( ai1, ai2+k1a*h/2.0, vi1, vi2+k1v*h/2.0, l1, l2, m1, m2);

  k3a = double2a( vi2+k2v*h/2.0);
  k3v = double2v( ai1, ai2+k2a*h/2.0, vi1, vi2+k2v*h/2.0, l1, l2, m1, m2 );

  k4a = double2a(vi2+k3v*h);
  k4v = double2v(ai1, ai2+k3a*h, vi1, vi2+k3v*h, l1, l2, m1, m2);

 
  *af2 = ai2 + (k1a + 2.0*(k2a+k3a) + k4a)* h/6.0;
  *vf2 = vi2 + (k1v + 2.0*(k2v+k3v) + k4v)* h/6.0;
}



int main(void)
{
  //Firstly I intend to explain to the user what this program does.
  printf("\n\n ----------------------------------------------\n");
  printf("|                                              |\n");
  printf("|      WE WILL BE COMPARING THE MOTION OF A    |\n");
  printf("|      SINGLE AND DOUBLE PENDULUM IN THIS      |\n");
  printf("|      PROGRAM TO DEMONSTRATE CHAOS IN PYSICS  |\n");
  printf("|                                              |\n");
  printf(" ----------------------------------------------\n\n\n"); 
 
  //First the single pendulum
  printf("----------------------------------------\n");
  printf("                    |\\ \n");
  printf("                    |A\\ \n");
  printf("                    |  \\ \n");
  printf("                        \\L\n");
  printf("                         \\ \n");
  printf("                          \\ \n");
  printf("                           \\ \n");
  printf("                           (M) \n\n");
  printf("Above is what a single pendulum looks like.\n");
  
  //Now we need to get inputs for the program

  //Get A
  printf("'A' represents angle of the pendulum in degrees.\n");
  printf("Please enter a  starting value of A (between 10 and 179):");

  int tempbool = 0;
  double angle = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&angle);
      if(angle < 10)
	{
	printf("\nThat is a small angle, choose a larger one for the motion to be more pronounced:");
	tempbool =0;
	}
	else if(angle >179)
	  {
	    printf("\nEnter an engle between 10 and 179 please:");
	    tempbool = 0;
	  }
    }
  //Convert angle to radians
  angle = angle*(PI/180);
  //Save to a temp in order to use later in double pendulm
  double tempa = angle;
  //Get L from the user
  printf("\n'L' represents the length of a rigid, string/rod in cm.\n");
  printf("Please enter the value of L:");
  tempbool = 0;
  double length = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&length);
      if(length <= 0)
	{
	printf("\nThe length can't be negative or 0. Try again:");
	tempbool =0;
	}
    }

  //Get M from the user
  printf("\n'M' represents the mass of the pendulum ball/string in grams.");
  printf("\nPlease enter the value of M:");
  tempbool = 0;
  double mass = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&mass);
      if(mass <= 0)
	{
	printf("\nThe mass can't be negative or 0. Try again:");
	tempbool =0;
	}
    }
  printf("\nOur angle in radians is: %lf\n Our length in cm is: %lf\n Our mass in g is: %lf\n",angle,length,mass);


  //The calling of the RK4 method for the single pendulum
  double ti=0.0; //Initialize t to 0
  double ai=angle; //Initialize x to the angle 
  double vi=0.0; //Initialize the velocity to 0
  
  double dt=0.01;    //step size
  double tmax = 5*2*PI*sqrt(length/GRAV); //Set tmax to 5 periods of the pendulum

  double tf, af, vf;
  tf=ti;
  af=ai;
  vf=vi;

  //Prepare to write to file:
  FILE* fPtr;
  char fileName[] = "spendulum.dat";
  fPtr = fopen(fileName, "w");

  fprintf(fPtr,"%lf  %lf  %lf  %lf \n",tf,X(af, length), Y(af, length), vf);

  while(ti<tmax)
    {
    tf += dt; //increment t

    rk4(ti, ai, vi, tf, &af, &vf, length);

    fprintf(fPtr,"%lf  %lf  %lf  %lf \n",tf,X(af, length), Y(af, length),vf);
    //now update the vars
    ti=tf; 
    ai=af;
    vi=vf;
    }
  fclose(fPtr);
  printf("\n\nTHE DATA HAS BEEN SAVED TO 'spendulum.dat'!");

  //Now lets do the double pendulum to demonstrate Chaos!
  printf("\n\n ----------------------------------------------\n");
  printf("|                                              |\n");
  printf("|      NOW WE WILL MODEL A DOUBLE PENDULUM     |\n");
  printf("|        AND COMPARE THE RESULTS TO ONE        |\n");
  printf("|                SINGULAR PENDULUM             |\n"); 
  printf("|                                              |\n");
  printf(" ----------------------------------------------\n\n\n"); 
 
  //Now show what a double pendulum looks like
  printf("----------------------------------------\n");
  printf("                    |\\ \n");
  printf("                    | \\ \n");
  printf("                    |A1\\ \n");
  printf("                    |   \\L1\n");
  printf("                         \\ \n");
  printf("                          \\ \n");
  printf("                           \\ \n");
  printf("                           (M1) \n");
  printf("                           /|   \n");
  printf("                          / |    \n");
  printf("                       L2/A2|     \n");
  printf("                        /   |      \n");
  printf("                       /            \n");
  printf("                      /              \n");
  printf("                     /                \n");
  printf("                   (M2)                 \n\n");
  printf("Above is what a double pendulum looks like.\n");
  
  //Reset angle to original
  angle = tempa;
  //Get A2 from the user
  printf("\n'A2' represents the angle of the second pendulum .\n");
  printf("Please enter the initial value of A2:");

  tempbool = 0;
  double angle2 = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&angle2);
      if(angle2 < 10 && angle2 > 10)
	{
	printf("\nThat is a small angle, choose a larger one for the motion to be more pronounced:");
	tempbool =0;
	}
	else if(angle2 >179 || angle2 < -179)
	  {
	    printf("\nEnter an engle between -179 and 179 please:");
	    tempbool = 0;
	  }
    }
  printf("\n'L2' represents the length of the second  rigid, string/rod in cm.\n");
  printf("Please enter the value of L2:");
  tempbool = 0;
  double length2 = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&length2);
      if(length <= 0)
	{
	printf("\nThe length can't be negative or 0. Try again:");
	tempbool =0;
	}
    }
  printf("\n'M2' represents the mass of the second pendulum ball/string in grams.");
  printf("\nPlease enter the value of M2:");
  tempbool = 0;
  double mass2 = 0;   
  while(tempbool == 0)
    {
      tempbool = 1;
      scanf("%lf",&mass2);
      if(mass2 <= 0)
	{
	printf("\nThe mass can't be negative or 0. Try again:");
	tempbool =0;
	}
    }
  printf("\nFor the second pendulum, our angle in radians is: %lf\n Our length in cm is: %lf\n Our mass in g is: %lf\n",angle2,length2,mass2);



  //The calling of the RK4 method for the double pendulum
  ti=0.0; //Initialize t to 0
  double ai1=angle; //Initialize x to the angle 
  double vi1=0.0; //Initialize the velocity to 0
  //need another set of vars for second pendulum
  double ai2 = angle2;
  double vi2 = 0.0; 
 

  double af1, vf1, af2, vf2;
  tf=ti;
  af1=ai1;
  vf1=vi1;
  af2 = ai2;
  vf2 = vi2;
  
  //Prepare to write to file:
  FILE* fPtr2;
  char fileName2[] = "dpendulum.dat";
  fPtr2 = fopen(fileName2, "w");

  fprintf(fPtr2,"%lf  %lf  %lf  %lf  %lf\n",tf,X(af1, length), Y(af1, length),X(af2, length2)+X(af1, length), Y(af2, length2)+Y(af1, length));

  while(ti<tmax)
    {
    tf += dt; //increment t

    doublerk4(ti, ai1, vi1, ai2, vi2, tf, 
	      &af1, &vf1, &af2, &vf2, length, length2, mass, mass2);


    fprintf(fPtr2,"%lf  %lf  %lf  %lf  %lf\n",tf,X(af1, length), Y(af1, length),X(af2, length2)+X(af1, length), Y(af2, length2)+Y(af1, length));
    //now update the vars
    ti=tf; 
    ai1=af1;
    vi1=vf1;
    ai2=af2;
    vi2 = vf2;
    }
  fclose(fPtr2);
  printf("\n\nTHE DATA HAS BEEN SAVED TO 'dpendulum.dat'!\n\n");
  return 0;
}
