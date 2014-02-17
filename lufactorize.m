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
#Exercise 2 matrices, LU factorized

a1 = [3, 1, 2; 6, 3, 4; 3, 1, 5];
[l1, u1] = lufactor(a1, 3);

disp("L:");
disp(l1);
disp("\nU:");
disp(u1);

disp("\nCheck: L * U");
disp(l1 * u1);

a2 = [4, 2, 0; 4, 4, 2; 2, 2, 3];
[l2, u2] = lufactor(a2, 3);

disp("L:");
disp(l2);
disp("\nU:");
disp(u2);

disp("\nCheck: L * U");
disp(l2 * u2);

a3 = [1, -1, 1, 2; 0, 2, 1, 0; 1, 3, 4, 4; 0, 2, 1, -1];
[l3, u3] = lufactor(a3, 4);

disp("L:");
disp(l3);
disp("\nU:");
disp(u3);

disp("\nCheck: L * U");
disp(l3 * u3);

#Exercise 4 systems, solved LU-style
disp("\nSystems:");
a1 = [3, 1, 2; 6, 3, 4; 3, 1, 5];

b1 = [0, 1, 3];
x1 = lusolve(a1, b1, 3);

disp(x1);
b2 = [2, 4, 6];
x2 = lusolve(a2, b2, 3);

disp(x2);

disp("2/4");
a4 = [1, 0, 0, 1; -1, 1, 0, 1; -1, -1, 1, 1; -1, -1, -1, 1];
[l4, u4] = lufactor(a4, 4);
disp(l4);
disp(u4);

disp("2/4 #5");
a5 = [1e-15, 1; 1, 1];
[l5, u5] = lufactor(a5, 2);
disp(l5);
disp("");
disp(u5);