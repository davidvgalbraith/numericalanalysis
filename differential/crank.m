% Program 8.4 Crank-Nicolson method
%with Dirichlet boundary conditions
% input: space interval [xl,xr], time interval [yb,yt],
%number of space steps M, number of time steps N
% output: solution w
% Example usage: w=crank(0,1,0,1,10,10)
function [x, t, w]=crank(xl,xr,yb,yt,M,N)
  f=@(x) e^(-x^2);
  l=@(t) 0*t;
  D=4.2*10^-5 * 3600;
  h=(xr-xl)/M;k=(yt-yb)/N;
  sigma=D*k/(h*h); m=M; n=N;
  a=diag(2+2*sigma*ones(m,1))+diag(-sigma*ones(m-1,1),1);
  a=sparse(a+diag(-sigma*ones(m-1,1),-1));
  a(m,:) = [zeros(1, m-3) -1 4 3];
  b=diag(2-2*sigma*ones(m,1))+diag(sigma*ones(m-1,1),1);
  b=sparse(b+diag(sigma*ones(m-1,1),-1));
  lside=ones(1, n+1);
  w(:,1)=etothe(xl+(1:m)*h)';
  for j=1:n
    sides=[lside(j)+lside(j+1);zeros(m-1,1)];
    w(:,j+1)=a\(b*w(:,j)+sigma*sides);
  end
  w=[lside;w];
  x=xl+(0:M)*h;t=yb+(0:N)*k;
  mesh(x,t,w');
  view (60,30); axis([xl xr yb yt -1 1])
end

function Y = etothe(X)
  for j=1:size(X)(1)
    for k=1:size(X)(2)
      Y(j,k) = e^(-X(j,k)^2);
    end
  end
end