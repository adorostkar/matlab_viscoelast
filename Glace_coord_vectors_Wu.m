% Initial quadrilateral mesh
% b.c. such that the solution is displ_x=0; displ_y - linear
%
% Face_flag(iface)=disco_region
%---->  Regular meshsize in each direction
% Data: Nx - number of intervals in x-direction
%       Ny - number of intervals in y-direction
%       [hx] - variable stepsize in x-direction
%       [hy] - variable stepsize in y-direction

function [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_Wu(L,H)
% ------- subdomain 1:

hx1 = L*1e-3;
xc1 =[0:hx1:L/5];  % till 2000 km
xc2 = [];
last= xc1(length(xc1));
q   = 1.02;
for k=1:41
    last = last*q;
    xc2  = [xc2 last];
end
last= xc2(41);
q   = 1.2;
xc3 = [];
for k=1:4
    last = last*q;
    xc3 = [xc3 last];
end
xc = [xc1 xc2 xc3 L];
Nx= length(xc);


hy1 = -1e-4;
yc1 =[0:hy1:-5*1e-4];
yc2 = [];
last= yc1(length(yc1));
q   = 1.5;
for k=1:2
    last = last*q;
    yc2  = [yc2 last];
end
yc2  = [yc2 -18e-4];
yc3 =[-18e-4+hy1:hy1:-100e-4];
last= -100e-4;
q   = 1.2; %1.06;
yc4 = [];
for k=1:12
    last = last*q;
    yc4 = [yc4 last];
end
yc4 = [yc4 -1000e-4];
last= -1000e-4;
q   = 1.5;
yc5 = [];
for k=1:5
    last = last*q;
    yc5 = [yc5 last];
end
yc = [yc1 yc2 yc3 yc4 yc5 H];
Ny= length(yc);

hx = 1e-3;
hy = 1e-4;
return
