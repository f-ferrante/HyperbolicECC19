
function [x1,t1, V1to, V2to]=hyp_static_bc
% Code originally written by Ying Tang to illustrate a CDC2014 paper
% and then adapted by Francesco Ferrante to illustrate a CDC2018 paper. The
% codes relis on the software for the solution to hyperbolic PDEs written
% by Shampine and freely available at http://faculty.smu.edu/shampine/current.html

global  K H lamb1 lamb2 Thorizon Lambda
lamb1 = Lambda(1,1);
lamb2 = Lambda(2,2);
Thorizon =30;
method = 'LxF';
pdefun = @pdes;

% Mesh grid
Nx = 20; %number of points in the space domain   
x = linspace(0,1,Nx);
Tinit = 0;
t = Tinit;
% Initial conditions %%%% SHOULD SATISFY THE COMPATIBILITY CONDITION
V(1,:)=(cos(4*pi*x)-1);
V(2,:)=(cos(2*pi*x)-1);


sol = setup(1,pdefun,t,x,V,method,[],@bcfun);

% check CFL condition
dx = x(2)-x(1);
dt = 0.95*dx/10;

Nt=floor(Thorizon/dt);
howfar = Thorizon/Nt;

%to stack the solution
V1to=zeros(Nt,Nx);
V2to =zeros(Nt,Nx);

V1to(1,:) = V(1,:);
V2to(1,:) = V(2,:);

for m = 1:Nt-1
    disp([num2str(m) '/' num2str(Nt)]) 
    %keyboard
    sol = hpde(sol,howfar,dt);  
    t = sol.t;
    V = sol.u;
    V1to(m+1,:) = V(1,:);
    V2to(m+1,:) = V(2,:);
   


end

%%%% Generate the two vectors z and t
t = linspace(Tinit,Thorizon,Nt);
[x1,t1] = meshgrid(x,t);
end

    
%=========================================================================
% Subfunctions

 function d = delta(x,t)
 global Thorizon;
if t<Thorizon/2
            d = 1+0*sin(x)*sin(1*t);
        else
            d=0;
end    
end
    % end function delta
    
    function F = pdes(t,x,V,V_x) % linear PDE
        global lamb1 lamb2
        F = zeros(size(V));
        F(1,:) = -lamb1.*V_x(1,:)+delta(x,t);
        F(2,:) = -lamb2.*V_x(2,:)+delta(x,t);
    % end function pdes
    end
   
    function [XL,XR] = bcfun(t,XLex,XRex) 
    %boundary conditions
    global K B H
    XLtmp = XLex;
    Xbc=(H+B*K)*[XRex(1);XRex(2)];
    XLtmp(1) = [1 0]*Xbc;
    XLtmp(2) = [0 1]*Xbc;
    XRtmp = XRex;
    XL=XLtmp;
    XR=XRtmp;
   end