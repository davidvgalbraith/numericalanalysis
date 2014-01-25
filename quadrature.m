function [ Int, fcnt, info ] = qad( fun, a, b, tol )
%Integrates fun from a to b to accuracy within tol.
%   by David Galbraith
h = (b - a) / 4;
x1 = fun(a);
x2 = fun(a + h);
x3 = fun(a + 2*h);
x4 = fun(a + 3*h);
x5 = fun(b);

fcnt = 5;
q1 = (b - a) / 6 * (x1 + 4 * x3 + x5);
q2 = h / 3 * (x1 + 2 * x3 + 4 * (x2 + x4) + x5);
q = q2 + 1/15 * (q2 - q1);
if (abs(q - q2) < tol)
    Int = q;
    info = 'Success!!';
    return
else
    ee = (a + b) / 2;
    [x, y] = sooperqad(x1, x2, x3, a, ee, fun, tol);
    [w, r] = sooperqad(x3, x4, x5, ee, b, fun, tol);
    Int = x + w;
    fcnt = fcnt + y + r;
    info = 'Success!!!';
    return
end

end

function [Int, fcnt] = sooperqad(x1, x3, x5, a, b, fun, tol)
% recursion.
h = (b - a) / 4;
x2 = fun(a + h);
x4 = fun(b - h);
fcnt = 2;
q1 = (b - a) / 6 * (x1 + 4 * x3 + x5);
q2 = h / 3 * (x1 + 2 * x3 + 4 * (x2 + x4) + x5);
q = q2 + 1/15 * (q2 - q1);
if (abs(q - q2) < tol)
    Int = q;
    return
else 
    ee = (a + b) / 2;
    [x, y] = sooperqad(x1, x2, x3, a, ee, fun, tol);
    [w, r] = sooperqad(x3, x4, x5, ee, b, fun, tol);
    Int = x + w;
    fcnt = fcnt + y + r;
    return
end
end
