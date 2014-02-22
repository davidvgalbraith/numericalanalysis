%Computes the QR Factorization two ways

%Gram shmit method

function [Q, R] = grams(A, m, n)
  for j = 0 : n-1
    y = A(j*m+1:j*m+m);
    Aj = A(j*m+1:j*m+m);
    for i = 0 : j-1
      R(i+1, j+1) = Q(i * m + 1: i * m + m).' * Aj;
      y = y - R(i + 1, j + 1) * Q(i * m + 1: i * m + m);
    end
    R(j + 1, j + 1) = 2norm(y, m);
    Q(j*m+1:j*m+m) = y / R(j + 1, j + 1);
  end
end

A = [1, -4; 2, 3; 2, 2];

[Q, R] = grams(A, 3, 2);

disp(Q);

disp(R);