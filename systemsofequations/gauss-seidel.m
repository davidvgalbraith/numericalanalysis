#infinity vector norm
function x = infvnorm(a, n)
  x = 0;
  for m = 1:n
    if abs(a(m)) > x
      x = abs(a(m));
    end
  end
end


function [x, graphx, graphy] = gausidel(A, b, n, tol)
  for a = 1:n-1
    D(a, a) = A(a, a);
    Dinv(a, a) = 1.0/A(a, a);
    for bb = a+1:n
      U(a, bb) = A(a, bb);
      L(bb, a) = A(bb, a);
    end
    U(a, a+1) = A(a, a+1);
    L(a+1, a) = A(a+1, a);
  end
  D(n, n) = A(n, n);
  U(n, n) = 0;
  L(n, n) = 0;
  Dinv(n, n) = 1.0/A(n, n);
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