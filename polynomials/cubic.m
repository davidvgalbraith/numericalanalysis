
function [ rts, info ] = cubic21852365( a, b, c, d )
% find the rooots
% abcd are the coefficients of the degree-three polynomial whose zeros this routine finds.
if (a > 1000000000)
    b = b / a;
    c = c / a;
    d = d / a;
    a = 1;
end
if (abs(a) > 0)
if (abs(a) < 1.0e-8)
    b = b / a;
    c = c / a;
    d = d / a;
    a = 1;
end
end
if a == 0
    if b == 0
        if c == 0
            if d == 0
                r1 = 6 / 0;
                r2 = 7 / 0;
                r3 = 8 / 0;
                rts = [r1 r2 r3];
                info = 'HELLA ROOOOTS';
return
        else
            info = 'NO ROOOOTS';
            rts = [0 / 0];
return
            end
        else
            r1 = -d / c;
            rts = [r1];
            info = 'ONE ROOOOT';
return
        end
    else
        if (abs(b) < 1e-8)
            c = c / b;
            d = d / b;
            b = 1;
        end
        r1 = posquart(b, c, d);
        r2 = negsquart(b, c, d);
        if (r1 == r2)
            rts = [r1];
            info = 'ONE ROOOOOT TWO MULTAPLICITY';
            return
        else
        info = 'TWO ROOOOTS';
        rts = [r1 r2];
return
        end
    end
else
r1 = newt(a, b, c, d);
if b == 0 && c == 0 && d == 0
    rts = [0];
    info = 'ONE ROOOOT';
    return
end
e = [1 -r1];
f = [a b c d];
[q, r] = deconv(f, e);
q1 = dot([1 0 0], q);
q2 = dot([0 1 0], q);
q3 = dot([0 0 1], q);
r2 = posquart(q1, q2, q3);
r3 = negsquart(q1, q2, q3);
for kk = 1:50
    r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
    r2 = r2 - eval(a, b, c, d, r2) / eval(0, 3 * a, 2 * b, c, r2);
    r3 = r3 - eval(a, b, c, d, r3) / eval(0, 3 * a, 2 * b, c, r3);
end
rts = [r1 r2 r3];
info = 'THREE ROOOOTS';
        
end
return

function r1 = posquart( a, b, c )
if (b < 0)
r1 = (-2 * c) / (b + sqrt(b * b - 4 * a * c));
else
    r1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
end
return

function r2 = negsquart( a, b, c )
if (b < 0)
r2 = (-2 * c) / (b - sqrt(b * b - 4 * a * c));
else
    r2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
end
return


function r1 = newt(a, b, c, d)
r1 = 987654321;
r = 1.0 / a * (abs(b) + abs(c) + abs(d));
if abs(r) < 1.0
r = 1.0;
end
e = -r;
if (r > 1000)
    e = -r;
    w = 0;
    while((r - e) > 4e-6)
        if (sign(eval(a, b, c, d, w)) == sign(eval(a, b, c, d, e)))
            e = w;
            w = (e + r) / 2;
        else
            r = w;
            w = (e + r) / 2;
        end
    end
end
if (e > r)
    e = -e;
    r = -r;
end
for s = e: .01 : r
t = s;
for n = 1:50
t = t - eval(a, b, c, d, t) / eval(0, 3 * a, 2 * b, c, t);
if(abs(eval(a, b, c, d, s)) < 1)
end
if abs(eval(a, b, c, d, t)) < eps * 50
r1 = t - eval(a, b, c, d, t) / eval(0, 3 * a, 2 * b, c, t);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
r1 = r1 - eval(a, b, c, d, r1) / eval(0, 3 * a, 2 * b, c, r1);
return
end
end
end
if (r1 == 987654321)
    r1 = (r + e) / 2;
    return
end



function y = eval(a, b, c, d, x)
y = d + x *(c + x * (b + x * a));
return
