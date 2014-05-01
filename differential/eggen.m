function s=eggen(Nx)
  for s=0:0.01:2
    %% Lax-Wendroff
    Alw = (1-s^2)*eye(Nx)...
        +(1+s)/2*s*diag(ones(Nx-1,1),-1)...
        +(s-1)/2*s*diag(ones(Nx-1,1), 1);

    %% Lax-Friedrichs
    Alf = zeros(Nx)...
        +(1+s)/2*diag(ones(Nx-1,1),-1)...
        +(1-s)/2*diag(ones(Nx-1,1), 1);

    %% Crank-Nicolson
    Acn1 = eye(Nx)...
        -0.25*s*diag(ones(Nx-1,1),-1)...
        +0.25*s*diag(ones(Nx-1,1), 1);
    Bcn = eye(Nx)...
        +0.25*s*diag(ones(Nx-1,1),-1)...
        -0.25*s*diag(ones(Nx-1,1), 1);
    Acn = Acn1\Bcn; %% note that A^{-1}B is the matrix   that determines stability of CN scheme.
    
    L1 = eig(Alw); %% substitute matrix from scheme
    [useless,ii] = max(abs(L1));
    eval1 = L1(ii);
    %hold on, plot(real(eval1),imag(eval1),'y.','MarkerSize',20)
    L2 = eig(Alf); %% substitute matrix from scheme
    [junk,iii] = max(abs(L2));
    eval2 = L2(iii);
    %hold on, plot(real(eval2),imag(eval2),'g.','MarkerSize',20)
    L3 = eig(Acn); %% substitute matrix from scheme
    [crap,iiii] = max(abs(L3));
    eval3 = L3(iiii);
    hold on, plot(real(eval3),imag(eval3),'b.','MarkerSize',20)
    x = -1:.01:1;
    y = sqrt(1-x.^2);
    %hold on, plot(x, y,'r','LineWidth',2)
    %hold on, plot(x,-y,'r','LineWidth',2)
    axis equal % this makes the circle look like a circle

  end
end