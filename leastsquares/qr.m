%Computes the QR Factorization two ways

%Two norm of the m-dimensional vector y
function norm = twonorm(y, m)
  sum = 0;
  for i = 1:m
    sum += y(i) * y(i);
  end
  norm = sqrt(sum);
end
  

%Gram shmit method

function [Q, R] = grams(A, m, n)
  R = zeros(n);
  Q = zeros(m, n);
  for j = 0 : n-1
    y = A(j*m+1:j*m+m);
    Aj = A(j*m+1:j*m+m); 
    for i = 0 : j-1
      qi =  Q(i * m + 1: i * m + m);
      R(i+1, j+1) = qi * Aj.';
      y = y - R(i + 1, j + 1) * Q(i * m + 1: i * m + m);
    end
    R(j + 1, j + 1) = twonorm(y, m);
    Q(j*m+1:j*m+m) = y / R(j + 1, j + 1);
  end
end

%Makes A, an nxn matrix, the lower-right submatrix of the mxm identity

function B = merge(A, n, m) 
  B = zeros(m);
  for j = 1:m-n
    B(j, j) = 1;
  end
  for k = m-n+1:m
    for l = m-n+1:m
      B(k, l) = A(k - m + n, l - m + n);
    end
  end
end

%Householder method
function [Q, R] = house(A, m, n)
  tmp = A;
  R = zeros(n);
  Q = zeros(m, n);
  H = eye(m);
  for j = 0 : n-1
    x = A(j*m+1+j:j*m+m);
    w = zeros(1, m-j);
    w(1) = twonorm(x, m-j);
    v = w - x;
    vv = (eye(m-j) - 2 * (v.' * v) / (v * v.'));
    Hn = merge(vv, m-j, m);
    H = Hn * H;
    A = H * A;
  end
  R = H * tmp;
  Q = inv(H);
end
