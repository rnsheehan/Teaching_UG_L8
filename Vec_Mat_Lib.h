#ifndef VEC_MAT_LIB_H
#define VEC_MAT_LIB_H

//#ifndef ATTACH_H
//#include "Attach.h"
//#endif

// This library contains functions that dynamically create vectors and matrices of type double, of arbitrary size
// Indexing of arrays is one-based, i.e. indices start at 1 and end at N
// This is different from the C/C++ default of zero-based array indexing, i.e. indices start at zero and end at N-1
// R. Sheehan 31 - 1 - 2013

namespace vec_mat_funcs{

	double *vector(int size); // this function will return a vector that can hold size elements of type double

	double inf_norm(double *x, int size); // compute the infinity norm of a vector

	double *vec_diff(double *a, double *b, int size); // compute the difference between two vectors

	double **matrix(int rows, int columns); // this function will return a matrix that can hold row*cols elements of type double

	void print_vector(double *vec, int size); // print a vector to the screen

	void print_matrix(double **mat, int rows, int columns); // print a matrix to the screen

}

#endif