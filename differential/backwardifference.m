				% Program 8.2 Backward difference method for heat equation
				% input: space interval [xl,xr], time interval [yb,yt],
				% number of space steps M, number of time steps N
				% output: solution w
				% Example usage: w=heatbd(0,1,0,1,10,10)

				%copyed from sauer
function w=heatbd(xl,xr,yb,yt,M,N)
  f=@(x) sin(pi*x);
  l=@(t) 0*t;
  r=@(t) 0*t;
  D=1/pi^2;
  h=(xr-xl)/M; k=(yt-yb)/N; m=M-1; n=N;
  sigma=D*k/(h*h);
  a=diag(1+2*sigma*ones(m,1))+diag(-sigma*ones(m-1,1),1);
  a=a+diag(-sigma*ones(m-1,1),-1);
  lside=l(yb+(0:n)*k); rside=r(yb+(0:n)*k);
  w(:,1)=f(xl+(1:m)*h)';
  for j=1:n
    w(:,j+1)=a\(w(:,j)+sigma*[lside(j);zeros(m-2,1);rside(j)]);
  end
  w=[lside;w;rside];
  x=(0:m+1)*h;t=(0:n)*k;
  mesh(x,t,w')
  view(60,30);axis([xl xr yb yt -1 2])
end