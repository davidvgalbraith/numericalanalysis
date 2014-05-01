function w=leapfrog(xl,xr,yb,yt,M,N)
  l=@(t) 0 * t;
  r = @(t) 0 * t;
  D=2;
  h=(xr-xl)/M; k=(yt-yb)/N; m=M-1; n=N;
  sigma=D*k/(h);
  a=diag(2-2*sigma*sigma * ones(m, 1))+diag(sigma*sigma * ones(m-1,1),1);
  a=a+diag(sigma*sigma*ones(m-1,1),-1);
  xes = xl + (1:m) * h;
  ts = yb + (0:n) * k;
  lside=l(ts); rside = r(ts);
  w(:, 1) = f(xes);
  rupert = 1/2 * a * f(xes)' + k * g(xes)';
  w(:, 2) = 1/2 * a * f(xes)' + k * g(xes)';
  for j=2:n
    w(:,j+1)=a*w(:,j)-w(:,j-1);
  end
  w = [lside; w; rside];
  x=(0:m+1)*h;t=(0:n)*k;
  mesh(x,t,w')
  view(60,30);axis([xl xr yb yt -1 1])
end

function y = f(xes)
  for a = 1:size(xes)(2)
    y(a) = 0;
  end
end

function z = g(xes)
  for a = 1:size(xes)(2)
    z(a) = 2 * pi * sin(pi * xes(a));
  end
end