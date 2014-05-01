#simultanius iteration
#evals is eigenvalues of A
#evecs is eigenvectors of A
#A is a square matrix
#k is number of iterations to proceed for
function [evals, evecs] = simulit(A, k)
  [m, n] = size(A);
  evecs = eye(m, m);
  for j=1:k
    [evecs, R] = qr(A*evecs);
  end
  evals = diag(evecs' * A * evecs);
end

disp("Simultaneous Iteration eigenvalues:");
A = [7, -33, -15; 2, 26, 7; -4, -50, -13];
[val, vec] = simulit(A, 1000);
val 
disp("Simultaneous Iteration eigenvectors");
vec
