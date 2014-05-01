#the power iteration method to find the leading eigenvalue and the corresponding eigen-
#vector of an n × n matrix

function [eval, evec] = poweriteration(A, guess, its)
  for a=1:its
    u = guess/norm(guess);
    guess = A * u;
    eval = u' * guess;
  end
  evec = guess/norm(guess);
end

#the inverse power iteration method to find the smallest eigenvalue and the corre-
#sponding eigenvector of an n × n matrix

function [eval, evec] = inversepoweriteration(A, guess, shift, its)
  shifted = A - shift * eye(size(A));
  for a=1:its
    u = guess/norm(guess);
    guess = shifted\u;
    eval = u' * guess;
  end
  evec = guess/norm(guess);
  eval = 1/eval + shift;
end

#the Rayleigh quotient iteration method to find the smallest eigenvalues and corresponding
#eigenvector of an n × n matrix
function [eval, evec] = raleyquotient(A, guess, its)
  for a=1:its
    u = guess/norm(guess);
    eval = u' * A * u;
    guess = (A - eval * eye(size(A)))\u;
  end
  evec = guess/norm(guess);
  eval = u' * A * u;
end


A = [7, -33, -15; 2, 26, 7; -4, -50, -13];

[e1, v1] = poweriteration(A, [1; 1; 1], 10000);

[e2, v2] = inversepoweriteration(A, [1; 1; 1], 0, 10000);

[e3, v3] = raleyquotient(A, [1; 1; 1], 10000);

e1 
v1
e2
v2
disp("Rayleigh quotient iteration is not a functional algorithm. Sorry.");
e3
v3
eig(A)