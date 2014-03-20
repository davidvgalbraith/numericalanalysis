# Naive Gausian Eliminimation

function x = eliminimate(a, b, n)

  for j = 1 : n-1
    if abs(a(j, j)) < eps; error('zero pivot encountered'); end
    for i = j+1 : n
      mult = a(i, j)/a(j, j);
      for k = j+1 : n
	a(i, k) = a(i, k) - mult*a(j, k);
      end
      b(i) = b(i) - mult*b(j);
    end
  end
  for i = n : -1 : 1
    for j = i+1 : n
      b(i) = b(i) - a(i, j) * x(j);
    end
    x(i) = b(i)/a(i, i);
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

#Infinity matrix norm
function x = infmnorm(A, n);
  x = 0;
  for k = 1:n
    y = 0;
    for m = 1:n
      y += abs(A(k, m));
    end
    if y > x
      x = y;
    end
  end
end

#find values of k for which the approximation in 2.3.1 has no correct
#significant figures
#prints the numbers 13-25 then crashes on a zero pivot
#conclusion: n >= 13 -> no significant figures correct
function x = nosigfigs(n)
  for k = 1:n
    [x, H] = flilberfy(k);
    #no significant figures right--approximation of 1 off by more than
    #1/2
    if infvnorm(x'-ones(k, 1), k) > 0.5
      disp(k);
    end
  end
end
