#include "util.hpp"

void info() {
	std::cout << std::endl;
}

void print() {
	std::cerr << std::endl;
}

void printSLE(const gsl_matrix* A, const gsl_vector* x,
              const gsl_vector* b, const int N) {
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			std::cerr << std::showpos << std::setprecision(10) << std::scientific
			          << gsl_matrix_get(A, i, j) << "\t";
		}
		std::cerr << std::showpos << std::setprecision(10) << std::scientific 
		          << "*\t" << gsl_vector_get(x, i)
		          << "=\t" << gsl_vector_get(b, i) << "\n";
	}
	std::cerr << std::endl;
}

void checkSLE(const gsl_matrix* A, const gsl_vector* x, 
              const gsl_vector* b, const int N) {
//	TODO - make complete logging to the file
//	print("Checking the solution of the SLE ... ");

	// Checking on condition number of matrix A
	int signum;
	gsl_permutation *perm = gsl_permutation_alloc(N);
	gsl_matrix *LU = gsl_matrix_alloc(N, N);
	gsl_matrix *invA = gsl_matrix_alloc(N, N);
	gsl_matrix_memcpy(LU, A);
	gsl_linalg_LU_decomp(LU, perm, &signum);
	gsl_linalg_LU_invert(LU, perm, invA);
	gsl_matrix_free(LU);
	gsl_permutation_free(perm);
	gsl_vector *row = gsl_vector_alloc(N);
	double normA = 0;
	double normInvA = 0;
	for(int i = 0; i < N; i++) {
		gsl_matrix_get_row(row, A, i);
		double dasum = gsl_blas_dasum(row);
		if (dasum > normA) normA = dasum;
		gsl_matrix_get_row(row, invA, i);
		dasum = gsl_blas_dasum(row);
		if (dasum > normInvA) normInvA = dasum;
	}
	double conditionNumber = normA * normInvA;
	if (conditionNumber > 1000)
		//print("Condition number of matrix of SLE is ", conditionNumber);
	gsl_vector_free(row);
	
	// Checking on Ax == b
	gsl_vector *tmp = gsl_vector_alloc(N);
	gsl_vector_memcpy(tmp, b);
	// tmp = A*x - b, i.e. error
	gsl_blas_dgemv(CblasNoTrans, 1, A, x, -1, tmp);
	for(int i = 0; i < N; i++) {
		if (fabs(gsl_vector_get(tmp, i) / gsl_vector_get(b, i)) > 1e-8) {
			if (gsl_vector_get(b, i) == 0) {
				if (fabs(gsl_vector_get(tmp, i)) > 1e-8)
					print("Ax =", gsl_vector_get(tmp, i), "at string", i,
					      ". But b = 0 here.");
			} else {
				print("( Ax - b ) / b =", 
				      gsl_vector_get(tmp, i) / gsl_vector_get(b, i),
				      "at string", i);
			}
		}
	}
	
	// Checking on inv(A)b == x
	gsl_vector_memcpy(tmp, x);
	// tmp = inv(A)*b - x, i.e. error
	gsl_blas_dgemv(CblasNoTrans, 1, invA, b, -1, tmp);
	for(int i = 0; i < N; i++) {
		if (fabs(gsl_vector_get(tmp, i) / gsl_vector_get(x, i)) > 1e-8)
			print("( inv(A)b - x ) / x =",
			      gsl_vector_get(tmp, i) / gsl_vector_get(x, i),
			      "at string ", i);
	}
	
	gsl_vector_free(tmp);
	gsl_matrix_free(invA);
	
//	TODO - make complete logging to the file	
//	print("Checking is done.");
}

double arctan(const double& a, const double& b) {
	double beta = atan(a / b);
	if (b == 0) {
		if (a > 0) beta = M_PI / 2;
		if (a < 0) beta = - M_PI / 2;
		if (a == 0) print("Error in atan: 0/0");
	}
	return beta;
}
