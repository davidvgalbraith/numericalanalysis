#The Conjugate Gradient Method, approximation-style


function [xk, graphx, graphy] = conjgrad(A, b, n, tol)
  xk = -1 * ones(n, 1);
  dk = b - A * xk;
  rk = b - A * xk;
  for k = 1:n

    q = 0;
    for m = 1:n
      if abs(rk(m)) > q
	q = abs(rk(m));
      end
    end
    graphx(k) = k;
    graphy(k) = q;
    if (q < tol)
      return
    else
      tmp = rk;
      ak = rk.' * rk / (dk.' * A * dk);
      xk  = xk + ak * dk;
      rk = rk - ak * A * dk;
      bk = rk.' * rk / (tmp.' * tmp);
      dk = rk + bk * dk;
    end
  end
end

#infinity vector norm
function x = infvnorm(a, n)
  x = 0;
  for m = 1:n
    if abs(a(m)) > x
      x = abs(a(m));
    end
  end
end

function [x, graphx, graphy] = jacobi(A, b, n, tol)
  D = zeros(n);
  Dinv = zeros(n);
  L = zeros(n);
  U = zeros(n);
  for a = 1:n
    D(a, a) = A(a, a);
    Dinv(a, a) = 1.0/A(a, a);
    for bb = a+1:n
      U(a, bb) = A(a, bb);
      L(bb, a) = A(bb, a);
    end
  end
  x = zeros(n, 1);
  graphx(1) = 0;
  graphy(1) = infvnorm(A * x - b, n);
  k = 2;
  while infvnorm(A * x - b, n) > tol
    x = Dinv*(b - (L + U)*x);
    graphx(k) = k-1;
    graphy(k) = infvnorm(A * x - b, n);
    k = k + 1;
  end
end


function [x, graphx, graphy] = gausidel(A, b, n, tol)
  D = zeros(n);
  Dinv = zeros(n);
  L = zeros(n);
  U = zeros(n);
  for a = 1:n
    D(a, a) = A(a, a);
    Dinv(a, a) = 1.0/A(a, a);
    for bb = a+1:n
      U(a, bb) = A(a, bb);
      L(bb, a) = A(bb, a);
    end
  end
  x = zeros(n, 1);
  graphx(1) = 0;
  graphy(1) = infvnorm(A * x - b, n);
  k = 2;
  while infvnorm(A * x - b, n) > tol
    next = zeros(n, 1);
    for l = 1:n
      next(l) = (Dinv * (b - U * x - L * next))(l);
    end
    graphx(k) = k-1;
    graphy(k) = infvnorm(A * next - b, n);
    k = k + 1;
    x = next;
  end
end
  
				% Program 2.1 Sparse matrix setup
				% Input: n = size of system
				% Outputs: sparse matrix a, r.h.s. b
function [a,b] = sparsesetup(n)
  e = ones(n,1); n2=n/2;
  a = spdiags([-e 3*e -e],-1:1,n,n);
	% Entries of a
  a(n2+1,n2) = -1; a(n2,n2+1) = -1;
				% Fix up 2 entries
  b=zeros(n,1);				% Entries of r.h.s. b
  b(1)=2;b(n)=2;b(2:n-1)=1;
end

[A, b] = sparsesetup(8);
[xk, graphx, graphy] = gausidel(A, b, 8, 1e-3);
disp(graphy);
%disp(xk);
#plot(graphx, graphy);
#title("Norm by iteration");
#xlabel("Iteration");
#ylabel("Norm of residual");