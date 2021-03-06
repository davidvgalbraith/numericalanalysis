function p = fire(m)
  if nargin < 1
    m = size(get(gcf,'colormap'),1); 
  end
  cmap_mat=[
	    [0;150;250]'/255;
	    1 1 1
	    1 130/255 0
	    [200;50;0]'/255;
	    ];

				%interpolate values
  xin=linspace(0,1,m)';
  xorg=linspace(0,1,size(cmap_mat,1));

  p(:,1)=interp1(xorg,cmap_mat(:,1),xin,'linear');
  p(:,2)=interp1(xorg,cmap_mat(:,2),xin,'linear');
  p(:,3)=interp1(xorg,cmap_mat(:,3),xin,'linear');
  find(p<=0)
end


dt = 1e-2; t=0:dt:1;
dx = 1e-2; x=0:dx:1;
[X,T] = meshgrid(x,t);
ue = exp(-T).*sin(pi*X);
surf(X,T,ue)
colormap nicecolormap
shading interp
view([0 90])
set(gcf,’Color’,’w’)
title(’Exact solution’)
