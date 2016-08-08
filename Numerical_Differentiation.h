#ifndef NUMERICAL_DIFFERENTIATION_H
#define NUMERICAL_DIFFERENTIATION_H

// function for implementing numerical differentiation via Richardson extrapolation

namespace numer_diff{

	double dfridr(double (*func)(double), double x, double h, double *err); 

}

#endif