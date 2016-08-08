#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of functions used to compute the solution of the linear BVP y''(x) = p(x) y'(x) + q(x) y(x) + r(x)
// on the domain a <= x <= b with Dirichlet boundary conditions y(a) = alpha, y(b) = beta
// R. Sheehan 9 - 3 - 2014

void bvp_solver::linear_bvp_solver(int n_internal, double start, double end, double alpha, double beta, double tolerance,\
					   double (*p)(double x), double (*q)(double x), double (*r)(double x), double *pos, double *sol)
{
	// function for computing the solution to the linear BVP y''(x) = p(x) y'(x) + q(x) y(x) + r(x)
	// on the domain a <= x <= b with Dirichlet boundary conditions y(a) = alpha, y(b) = beta
	// The solution is computed at the internal nodes of the mesh
	// The Gauss-Seidel iteration method is used to compute the solution to within the specified tolerance
	// The mesh is stored in pos and the solution is stored in sol
	// R. Sheehan 9 - 3 - 2014

	// 1. compute the node spacing for the mesh
	double h = (end - start) / (static_cast<double>( n_internal + 1 )); 
	double hsqr = DSQR(h); 
	double halfh = 0.5*h; 
	double x_val; 

	cout<<"Step-size = "<<h<<endl<<endl;

	int n_total = n_internal + 2; 

	// 2. store the start and end points of the solution
	pos[1] = start; pos[n_total] = end; 
	sol[1] = alpha; sol[n_total] = beta; 

	// 3. Create a matrix and vectors for holding the system of equations 
	double **A = matrix(n_internal, n_internal); 

	double *d = vector(n_internal); 

	// 4. Construct the system of equations whose solution represents that of the BVP
	// We know the matrix is tri-diagonal so we assign the values to the matrix accordingly
	
	for(int i=1; i<=n_internal; i++){

		x_val = start + (i*h); // note this is different from what was used for the solution of IVP

		pos[i+1] = x_val; 

		for(int j=1; j<=n_internal; j++){
			
			if(i == j){ // Main Diagonal
				A[i][j] = 2.0 + hsqr*(*q)(x_val); 
			}
			else if(i == j-1 ){ // Super Diagonal
				A[i][j] = -1.0  + halfh*(*p)(x_val);
			}
			else if(i == j+1){ // Sub Diagonal 
				A[i][j] = -1.0  - halfh*(*p)(x_val); 
			}
			else{
				A[i][j] = 0.0; 
			}

		}
	}

	/*if(n_internal < 25){
		cout<<"The matrix A is \n";
		print_matrix(A, n_internal, n_internal); 
	}*/

	// We also need to fill the RHS vector
	for(int i=1; i<=n_internal; i++){
		
		x_val = start + (i*h); // note this is different from what was used for the solution of IVP

		if(i==1){
			// include left boundary condition
			d[i] = -hsqr*(*r)(x_val) + (1.0 + halfh*(*p)(x_val))*alpha; 
		}
		else if(i==n_internal){
			// include right boundary condition
			d[i] = -hsqr*(*r)(x_val) + (1.0 - halfh*(*p)(x_val))*beta; 
		}
		else{
			// points in between
			d[i] = -hsqr*(*r)(x_val); 
		}

	}

	/*if(n_internal < 25){
		cout<<"The RHS vector d is \n"; 
		print_vector(d, n_internal); 
	}*/

	// Construct the vector that will contain the computed solution
	double *internal_sol = vector(n_internal); 

	// zero the elements in the vector
	// it is also possible to initialise this vector with an approximate solution
	for(int i=1; i<=n_internal; i++){
		internal_sol[i] = 0.0; 
	}

	// Compute the solution of A . x = d using an interative solution method Gauss-Seidel iteration
	int max_iterations = 10000; 

	bool converged; 

	double error; 

	gauss_seidel_solve(A, internal_sol, d, n_internal, max_iterations, converged, tolerance, error); 

	/*cout<<"The solution of the system A . x = d is \n"; 
	print_vector(internal_sol, n_internal); */

	// Store the computed solution in the vector sol
	for(int i=1l; i<=n_internal; i++){
		sol[i+1] = internal_sol[i]; 
	}

	// Delete vectors and matrices when they are no longer required
	delete[] A; 
	delete[] d; 
	delete[] internal_sol; 
}

void bvp_solver::gauss_seidel_solve(double **a, double *x, double *b, int n, int max_iter, bool &solved, double tol, double &error)
{
	// Solve the system Ax=b by Gauss-Seidel iteration method
	// This has a faster convergence rate than the last technique
	// Unless you can determine omega for a specific system use GS
	// R. Sheehan 3 - 8 - 2011
	
	bool cgt=false;

	double *oldx = vector(n); 
	double *diffx = vector(n); 
		
	int n_iter = 1;
	double bi, mi, err;
		
	while(n_iter<max_iter){

		// store the initial approximation to the solution
		for(int i=1; i<=n; i++){
			oldx[i] = x[i]; 
		}

		// iteratively compute the solution vector x
		for(int i=1; i<=n; i++){

			bi=b[i]; 
			mi=a[i][i];

			for(int j=1; j<i; j++){

				bi -= a[i][j]*x[j];

			}

			for(int j=i+1; j<=n; j++){

				bi -= a[i][j]*oldx[j];

			}

			x[i] = bi/mi;
		}
			
		// test for convergence
		// compute x - x_{old} difference between the two is the error
		diffx = vec_diff(x, oldx, n); 

		// error is measured as the length of the difference vector ||x - x_{old}||_{\infty}
		err = inf_norm(diffx, n); 
			
		if(abs(err)<tol){ // solution has converged, stop iterating

			cout<<"\nGauss-Seidel Iteration Complete\nSolution converged in "<<n_iter<<" iterations\n";

			cout<<"Error = "<<abs(err)<<endl<<endl;

			error=abs(err);

			cgt=solved=true;

			break;
		}
			
		n_iter++; // solution has not converged, keep iterating
	}
		
	if(!cgt){ // solution did not converge

		cout<<"\nError: Gauss-Seidel Iteration\n";

		cout<<"Error: Solution did not converge in "<<max_iter<<" iterations\n";

		cout<<"Error = "<<abs(err)<<endl;

		error=abs(err);

		solved=false;
	}

	delete[] oldx; 
	delete[] diffx; 
}

void bvp_solver::linear_bvp_solver_TDMA(int n_internal, double start, double end, double alpha, double beta, double tolerance,\
					   double (*p)(double x), double (*q)(double x), double (*r)(double x), double *pos, double *sol)
{
	// function for computing the solution to the linear BVP y''(x) = p(x) y'(x) + q(x) y(x) + r(x)
	// on the domain a <= x <= b with Dirichlet boundary conditions y(a) = alpha, y(b) = beta
	// The solution is computed at the internal nodes of the mesh
	// The tri-diagonal matrix algorithm is used to compute the solution to within the specified tolerance
	// The mesh is stored in pos and the solution is stored in sol
	// R. Sheehan 1 - 4 - 2014

	// 1. compute the node spacing for the mesh
	double h = (end - start) / (static_cast<double>( n_internal + 1 )); 
	double hsqr = DSQR(h); 
	double halfh = 0.5*h; 
	double x_val; 

	cout<<"Step-size = "<<h<<endl<<endl;

	int n_total = n_internal + 2; 

	// 2. store the start and end points of the solution
	pos[1] = start; pos[n_total] = end; 
	sol[1] = alpha; sol[n_total] = beta; 

	// 3. Create a matrix and vectors for holding the system of equations 
	
	double *a = vector(n_internal); // sub-diagonal 
	double *b = vector(n_internal); // main diagonal 
	double *c = vector(n_internal); // super diagonal
	double *d = vector(n_internal); // rhs vector

	// 4. Construct the system of equations whose solution represents that of the BVP
	// We know the matrix is tri-diagonal so we assign the values to the matrix accordingly
	
	for(int i=1; i<=n_internal; i++){

		x_val = start + (i*h); // note this is different from what was used for the solution of IVP

		pos[i+1] = x_val; 

		if(i>1){ // sub-diagonal 
			a[i] = -1.0 - halfh*(*p)(x_val);
		}

		b[i] = 2.0 + hsqr*(*q)(x_val); // main diagonal

		if(i<n_internal){
			c[i] = -1.0 + halfh*(*p)(x_val); // super-diagonal
		}
	}

	// We also need to fill the RHS vector
	for(int i=1; i<=n_internal; i++){
		
		x_val = start + (i*h); // note this is different from what was used for the solution of IVP

		if(i==1){
			// include left boundary condition
			d[i] = -hsqr*(*r)(x_val) + (1.0 + halfh*(*p)(x_val))*alpha; 
		}
		else if(i==n_internal){
			// include right boundary condition
			d[i] = -hsqr*(*r)(x_val) + (1.0 - halfh*(*p)(x_val))*beta; 
		}
		else{
			// points in between
			d[i] = -hsqr*(*r)(x_val); 
		}

	}

	/*if(n_internal < 25){
		cout<<"The RHS vector d is \n"; 
		print_vector(d, n_internal); 
	}*/

	// Construct the vector that will contain the computed solution
	double *internal_sol = vector(n_internal); 

	// zero the elements in the vector
	// it is also possible to initialise this vector with an approximate solution
	for(int i=1; i<=n_internal; i++){
		internal_sol[i] = 0.0; 
	}

	// Compute the solution of A . x = d using a tri-diagonal matrix algorithm
	TDMA(a, b, c, d, internal_sol, n_internal); 

	/*cout<<"The solution of the system A . x = d is \n"; 
	print_vector(internal_sol, n_internal); */

	// Store the computed solution in the vector sol
	for(int i=1l; i<=n_internal; i++){
		sol[i+1] = internal_sol[i]; 
	}

	// Delete vectors and matrices when they are no longer required
	
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d; 
	delete[] internal_sol; 
}

void bvp_solver::TDMA(double *a, double *b, double *c, double *d, double *x, int size)
{
	// Compute the solution of a tri-diagonal system of equations
	// a = sub-diagonal, b = main-diagonal, c = super-diagonal, d = rhs vector
	// solution is stored in x
	// R. Sheehan 21 - 3 - 2014
	
	double tmp; 
	double *cpr = vector(size); 
	double *dpr = vector(size); 

	//cpr[size] = 0; // not necessary

	// 1. Make the transformations of c to c' and d to d'
	for(int i=1; i<=size; i++){
		
		if(i==1){
			cpr[1] = c[1] / b[1]; 
			
			dpr[1] = d[1] / b[1]; 
		}
		/*else if(i == size){

			tmp = b[i] - ( cpr[i-1] * a[i] );
			
			dpr[i] = ( d[i] - ( dpr[i-1] * a[i] ) ) / tmp; 

		}*/
		else{
			tmp = b[i] - ( cpr[i-1] * a[i] );

			cpr[i] = c[i] / tmp; // this sets cpr[size] = 0 when i = size
			
			dpr[i] = ( d[i] - ( dpr[i-1] * a[i] ) ) / tmp; 
		}
	}

	/*cout<<"The vector c' is \n"; 
	print_vector(cpr, size); 

	cout<<"The vector d' is \n"; 
	print_vector(dpr, size); */

	// 2. Compute the solution using back substitution
	int stp = size; 
	int lasstp; 
	while(stp > 0){
		
		if(stp == size){
			x[stp] = dpr[stp]; 
		}
		else{
			x[stp] = dpr[stp] - ( cpr[stp] * x[lasstp] );
		}

		lasstp = stp; 

		stp--; 
	}

	delete[] cpr; 
	delete[] dpr; 
}