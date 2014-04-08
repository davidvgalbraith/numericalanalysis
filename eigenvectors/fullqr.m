#Full QR method
#Finds eigenvalues of any matrix ever

#start with Upper Heisenberg Form
#A will be upper heisenberg
#H will be the matrix of the household reflectors

function [eval, evec] = inversepoweriteration(A, guess, shift, its)
  randy = rand()/10;
  shifted = A - (shift - randy) * eye(size(A));
  for a=1:20
    u = guess/norm(guess);
    guess = shifted\u;
    eval = u' * guess;
  end
  evec = guess/norm(guess);
  eval = 1/eval + shift - randy;
end

function [A, H] = heisen(A)
  [m, n] = size(A);
  H = zeros(m, m);
  for k=1 : m-2
    x = A(k+1:m, k);
    H(1:m-k, k) = -sign(x(1)+eps) * norm(x) * eye(m-k, 1) - x;
    H(1:m-k,k)=H(1:m-k,k)/norm(H(1:m-k,k));
    A(k+1:m,k:m)=A(k+1:m,k:m)-2*H(1:m-k,k)*H(1:m-k,k)'*A(k+1:m,k:m);
    A(1:m,k+1:m)=A(1:m,k+1:m)-2*A(:,k+1:m)*H(1:m-k,k)*H(1:m-k,k)';
  end
end

function evals = shiftedqr(A)
  m = size(A, 1);
  tol = 1e-14;
  maxcount = 500;
  evals = zeros(m, 1);
  n=m;
  while n>1
    ccount = 0;
    while max(abs(A(n, 1:n-1))) > tol & ccount < maxcount
      ccount = ccount + 1;
      mu = A(n, n);
      [q, r] = qr(A - mu*eye(n));
      A = r * q + mu * eye(n);
    end
    if ccount < maxcount
      evals(n) = A(n, n);
      n = n-1;
      A = A(1:n, 1:n);
    else
      disc=(A(n-1,n-1)-A(n,n))^2+4*A(n,n-1)*A(n-1,n);
      evals(n)=(A(n-1,n-1)+A(n,n)+sqrt(disc))/2;
      evals(n-1)=(A(n-1,n-1)+A(n,n)-sqrt(disc))/2;
      n=n-2;
      A=A(1:n,1:n);
    end
  end
  if n>0
    evals(1) = A(1,1);
  end
end

function val = fullqr(A)
  A = heisen(A);
  val = shiftedqr(A);
end

#SVD of your three by three A
function [U, S, Vt] = svd(A)
  B = [zeros(size(A)), A'; A, zeros(size(A))];
  val = fullqr(B, 1000);
  [val1, vec1] = inversepoweriteration(B, [1; 1; 1; 1; 1; 1], val(2), 10000);
  [val2, vec2] = inversepoweriteration(B, [1; 1; 1; 1; 1; 1], val(4), 10000);
  [val3, vec3] = inversepoweriteration(B, [1; 1; 1; 1; 1; 1], val(6), 10000);
  vone = vec1(1:3)/norm(vec1(1:3));
  vtwo = vec2(1:3)/norm(vec2(1:3));
  vthree = vec3(1:3)/norm(vec3(1:3));
  uone = A * vone / val1;
  utwo = A * vtwo / val2;
  uthree = A * vthree / val3;
  U = [uone, utwo, uthree];
  S = [val1, 0, 0; 0, val2, 0; 0, 0, val3];
  Vt = [vone, vtwo, vthree]';
end
disp("My eigenvalues, vectors of A");
A = [7, -33, -15; 2, 26, 7; -4, -50, -13];

val = fullqr(A, 1000);

[val1, vec1] = inversepoweriteration(A, [1; 1; 1], val(1), 10000);
[val2, vec2] = inversepoweriteration(A, [1; 1; 1], val(2), 10000);
[val3, vec3] = inversepoweriteration(A, [1; 1; 1], val(3), 10000);

val1
vec1
val2
vec2
val3
vec3

disp("My SVD decomposition of A")
[U, S, Vt] = svd(A)

disp("Check: U*S*Vt")

U * S * Vt