Numerical Analysis routines

systemsofequations:

conjgrad.m implements the Conjugate Gradient iterative method to find the solutions of a symmetric, positive-definite matrix system of equations within a given tolerance.

jacobi.m implements the Jacobi iterative method to find the solutions of a symmetric, positive-definite matrix system of equations within a given tolerance.

gauss-seidel.m implements the Gauss-Seidel iterative method to find the solutions of a symmetric, positive-definite matrix system of equations within a given tolerance.

naivegauss.m implements naive Gaussian elimination to solve a system of linear equations. 

lufactorize.m implements LU factorization, turning a matrix A into the product of matrices LU where L is lower-triangular with ones on the diagonal and U is upper-triangular.


eigenvalues:

poweriteration.m implements the Power Iteration method to find the maximum eigenvalue and its corresponding eigenvector of a matrix, given the matrix, an initial guess, and the number of steps to iterate for, as well as the Inverse Power Iteration Method to find the an arbitrary eigenvalue and its corresponding eigenvector of a matrix, given the matrix, an initial guess of the vector, an initial guess of the eigenvalue (known as the "shift"), and a number of iterations.

polynomials:

cubic.m: Find the roots of a third-degree polynomial.


integrals:

quadrature.m: Integrate a given function between two endpoints, with a specified accuracy.
