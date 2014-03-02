# An application of Newton-Gauss minimization to computing the weather

%Two norm of the 3-dimensional vector y
function norm = twonorm(y)
  sum = 0;
  for i = 1:3
    sum += y(i) * y(i);
  end
  norm = sqrt(sum);
end
  
function next = m(v)
  t = 0.05;
  next(1) = v(1) + 10 * (v(2)- v(1)) * t;
  next(2) = v(2)+ v(1) * (28 - v(3)) * t;
  next(3) = v(3) + (v(1) * v(2)- 8.0 / 3.0 * v(3)) * t;
end

xzero = [-3.48, -3.3, 22.7]; xone = m(xzero); xtwo = m(xone);
xzerotrue = [-3.6, -3.2, 22.1]; xonetrue = m(xzerotrue); xtwotrue = m(xonetrue);
muzero = [-3.1, -3.7, 24]; muone = m(muzero); mutwo = m(muone);

ezero = twonorm(xzero - xzerotrue) / twonorm(xzerotrue)
eone = twonorm(xone - xonetrue) / twonorm(xonetrue)
etwo = twonorm(xtwo - xtwotrue) / twonorm(xtwotrue)

ezerohat = twonorm(muzero - xzerotrue) / twonorm(xzerotrue)
eonehat = twonorm(muone - xonetrue) / twonorm(xonetrue)
etwohat = twonorm(mutwo - xtwotrue) / twonorm(xtwotrue)

%--output--%
ezero =  0.027411; eone =  0.025396; etwo =  0.024105
ezerohat =  0.089630; eonehat =  0.080376; etwohat =  0.076061