#The Conjugate Gradient Method, approximation-style

#infinity vector norm
function x = infvnorm(a, n)
  x = 0;
  for m = 1:n
    if abs(a(m)) > x
      x = abs(a(m));
    end
  end
end

function [xk, graphx, graphy] = conjgrad(A, b, n, tol)
  xk = ones(n, 1)
  dk = b - A * xk
  rk = b - A * xk
  for k = 1:n
    q = infvnorm(rk)
    graphx(k) = k
    graphy(k) = q
    if (q < tol)
      return
    else
      tmp = rk
      ak = rk.' * r / (dk.' * A * dk)
      xk  = xk + ak * dk
      rk = rk - ak * A * dk
      bk = rk.' * rk / (tmp.' * tmp)
      dk = rk + bk * dk
    end
  end
end


				% Program 2.1 Sparse matrix setup
				% Input: n = size of system
				% Outputs: sparse matrix a, r.h.s. b
function [a,b] = sparsesetup(n)
  e = ones(n,1); n2=n/2;
  a = spdiags([-e 3*e -e],-1:1,n,n);
				% Entries of a
  c=spdiags([e/2],0,n,n);c=fliplr(c);a=a+c;
  a(n2+1,n2) = -1; a(n2,n2+1) = -1;
				% Fix up 2 entries
  b=zeros(n,1);
				% Entries of r.h.s. b
  b(1)=2.5;b(n)=2.5;b(2:n-1)=1.5;b(n2:n2+1)=1;
end

[A, b] = sparsesetup(100000)
[xk, graphx, graphy] = conjgrad(A, b, 100000, 1e-6)
%plot(graphx, graphy)
title("Norm by iteration")
xlabel("Iteration")
ylabel("Norm of residual")
