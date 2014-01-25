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

#Exercise 2.1.1

a1 = [2, -2, -1; 4, 1, -2; -2, 1, -1];
b1 = [-2, 1, 3];

a2 = [1, 2, -1; 0, 3, 1; 2, -1, 1];
b2 = [2, 4, 2];

a3 = [2, 1, -4; 1, -1, 1; -1, 3, -2];
b3 = [-7, -2, 6];

n = 3;

x1 = eliminimate(a1, b1, n);
x2 = eliminimate(a2, b2, n);
x3 = eliminimate(a3, b3, n);

disp(x1);
disp(x2);
disp(x3);

#Exercise 2.1.2

function x = hilberfy(n)
  for i = 1:n
    for j = 1:n
      H(i, j) = 1/(i + j - 1);
    end
  end
  b = ones(n);
  x = eliminimate(H, b, n);
  return
end

x4 = hilberfy(2);
x5 = hilberfy(5);
x6 = hilberfy(10);

disp(x4);
disp(x5);
disp(x6);

#Exercise 2.3.1

#Like hilbert matrix but a little different
function [xc, H] = flilberfy(n)
  for i = 1:n
    for j = 1:n
      H(i, j) = 5/(i + 2 * j - 1);
    end
  end
  x = ones(n);
  b = H * x;
  xc = eliminimate(H, b, n);
  return
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


disp("\n 2.3.1:");
disp("\tn = 6:");
[x1, H6] = flilberfy(6);
disp("Approximate x");
disp(x1);
disp("Forward error infinity norm:");

r = H6 * (ones(6, 1) - x1');
b = H6 * ones(6, 1);
forward = infvnorm(ones(6, 1) - x1', 6);
disp(forward);
disp("Magnification factor:");
back = infvnorm(r, 6);
disp(forward/(back/infvnorm(b, 6)));
disp("Condition Number:");
disp(infmnorm(H6, 6) * infmnorm(inv(H6), 6));


disp("\n\tn= 10:");
[x2, H10] = flilberfy(10);
disp("Approximate x");
disp(x2);
disp("Forward error infinity norm:");

r = H10 * (ones(10, 1) - x2');
b = H10 * ones(10, 1);
forward = infvnorm(ones(10, 1) - x2', 10);
disp(forward);
disp("Magnification factor:");
back = infvnorm(r, 10);
disp(forward/(back/infvnorm(b, 10)));
disp("Condition Number:");
disp(infmnorm(H10, 10) * infmnorm(inv(H10), 10));

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

disp("\n2.3.5");
nosigfigs(100);