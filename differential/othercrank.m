				% Program 8.4 Crank-Nicolson method
				%with Dirichlet boundary conditions
				% input: space interval [xl,xr], time interval [yb,yt],
				%number of space steps M, number of time steps N
				% output: solution w
				% Example usage: w=crank(0,1,0,1,10,10)
function [w, exact] = othercrank(xl,xr,yb,yt,M,N)
  f=@(x) 1/sqrt(2 * pi) * e^(-x^2/2);
  l=@(t) 0*t;
  r=@(t) 0*t;
  D=1;
  h=(xr-xl)/M;k=(yt-yb)/N;
  sigma=D*k/(h*h); m=M-1; n=N;
  a=diag(2+2*sigma*ones(m,1))+diag(-sigma*ones(m-1,1),1);
  a=a+diag(-sigma*ones(m-1,1),-1);
  b=diag(2-2*sigma*ones(m,1))+diag(sigma*ones(m-1,1),1);
  b=b+diag(sigma*ones(m-1,1),-1);
  lside=l(yb+(0:n)*k); rside=r(yb+(0:n)*k);
  w(:,1)= 1/sqrt(2 * pi) * etothe(xl+(1:m)*h)';
  for j=1:n
    sides=[lside(j)+lside(j+1);zeros(m-2,1);rside(j)+rside(j+1)];
    w(:,j+1)=a\(b*w(:,j)+sigma*sides);
  end
  w=[lside;w;rside];
  x=xl+(0:M)*h;t=yb+(0:N)*k;
  [X, T] = meshgrid(x, t);
  exact = 1/sqrt(2 * pi) * etothe(X - T);
  disp(1);
  %mesh(x,t,w');
  disp(2);
  mesh(x, t, exact);
  disp(3);
  view (60,30); axis([xl xr yb yt -1 1])
end

function Y = etothe(X)
  for j=1:size(X)(1)
    for k=1:size(X)(2)
      Y(j,k) = e^(-X(j,k)^2);
    end
  end
end