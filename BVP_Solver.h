#ifndef BVP_SOLVER_H
#define BVP_SOLVER_H

// Declaration of functions used to compute the solution of the linear BVP y''(x) = p(x) y'(x) + q(x) y(x) + r(x)
// on the domain a <= x <= b with Dirichlet boundary conditions y(a) = alpha, y(b) = beta
// R. Sheehan 9 - 3 - 2014

namespace bvp_solver{

	// linear BVP solver
	void linear_bvp_solver(int n_internal, double start, double end, double alpha, double beta, double tolerance,\
						   double (*p)(double x), double (*q)(double x), double (*r)(double x), double *pos, double *sol); 

	// Gauss-Seidel iterative solver for the system A x = b
	void gauss_seidel_solve(double **a, double *x, double *b, int n, int max_iter, bool &solved, double tol, double &error); 

	// linear BVP solver with tri-diagonal matrix algorithm
	void linear_bvp_solver_TDMA(int n_internal, double start, double end, double alpha, double beta, double tolerance,\
						   double (*p)(double x), double (*q)(double x), double (*r)(double x), double *pos, double *sol);

	// Thomas algorithm for the solution of tri-diagonal systems
	void TDMA(double *a, double *b, double *c, double *d, double *x, int size); 

}

#endif