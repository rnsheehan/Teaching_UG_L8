#ifndef ATTACH_H
#include "Attach.h"
#endif

double numer_diff::dfridr(double (*func)(double),double x,double h,double *err)
{
	// Returns the derivative of a function func at the point x by Ridder's Method of polynomial extrapolation
	// Method is based on Richardson extrapolation
	// An estimate of the error is returned in err
	// Implementation taken from "Numerical Recipes in C" by Press et al. 

	int i,j;
	double errt,fac,hh,ans;
	double **a;

	//Constants necessary for the calculation
	static const double CON = 1.4;
	static const double CON2 = DSQR(CON);
	static const double SAFE = 2.0;
	static const double BIG = (1.0e30);
	static const int NTAB = 10;

	if(h==0.0){
		cerr<<"H must be non-zero in dfridr\n";	
	}

	a=matrix(NTAB,NTAB);

	hh=h;

	a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh); // this is the standard central difference approximation

	*err=BIG;

	for(i=2;i<=NTAB;i++){

		hh/=CON;

		a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh); // central difference approximation with smaller stepsize

		fac=CON2;

		for(j=2;j<=i;j++){// Compute extrapolations of various orders			

			a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1.0);

			fac=CON2*fac;

			errt = max(fabs(a[j][i] - a[j-1][i]), fabs(a[j][i] - a[j-1][i-1])); 

			if(errt <= *err){

				*err = errt; // If error is decreased then save the improved answer
				
				ans = a[j][i]; 
			}

		}

		if(fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err) ){

			delete[] a; // If higher order is worse by factor SAFE then quit

			return ans;
		}

	}

	delete[] a;
	return ans;
}