#compute the LU Factorization (2.2.1)

function [L, A] = lufactor(A, n)
   for j = 1 : n-1
    if abs(A(j, j)) < eps; error('zero pivot encountered'); end
    L(j, j) = 1;
    for i = j+1 : n
      mult = A(i, j)/A(j, j);
      L(i, j) = mult;
      for k = 1 : n
	A(i, k) = A(i, k) - mult*A(j, k);
      end
    end
  end
  L(n, n) = 1;
end

#Add back substitution to solve systems (2.2.2)

function x = lusolve(A, b, n)
  [L, U] = lufactor(A, n);
  c = frontsub(L, b, n);
  x = backsub(U, c, n);
  return
end

#back-substitution for Ax = b, A upper-triangular
function x = backsub(A, b, n)
  for i = n : -1 : 1
    for j = i+1 : n
      b(i) = b(i) - A(i, j) * x(j);
    end
    x(i) = b(i) / A(i, i);
  end
end

#back-substitution for Ax = b, A lower-triangular
function x = frontsub(A, b, n) 
  for i = 1 : n
    for j = 1 : i-1
      b(i) = b(i) - A(i, j) * x(j);
    end
    x(i) = b(i) / A(i, i);
  end
end
