% Initial quadrilateral mesh
% b.c. such that the solution is displ_x=0; displ_y - linear
%
% Face_flag(iface)=disco_region
%---->  Regular meshsize in each direction
% Data: Nx - number of intervals in x-direction
%       Ny - number of intervals in y-direction
%       [hx] - variable stepsize in x-direction
%       [hy] - variable stepsize in y-direction

function [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_Wu_coarse(L,Hl,l)
% ------- subdomain 1:
% x-coordinates
hx1 = l/10;
xc1 =[0:hx1:2*l];  % twice as the length of the ice
Lr = L-2*l;
Lrr=roundn(Lr/3,-2);
Lrl=Lr-2*Lrr;
hx2 = 2*hx1;
xc2=[xc1(end)+hx2:hx2:Lrr+l];
hx3 = 2*hx2;
xc3=[xc2(end)+hx3:hx3:2*Lrr+l];
hx4 = (L-xc3(end))/4;
xc4=[xc3(end)+hx4:hx4:L];

xc = [xc1 xc2 xc3 xc4];
Nx= length(xc);

% y-coordinates
H = Hl(end);
hy1 = -l/10;
yc1 =[0:hy1:-l];
Hr = H+l;
Hrr=roundn(Hr/3,-2);
Hrl=Hr-2*Hrr;
hy2 = 2*hy1;
yc2=[yc1(end)+hy2:hy2:Hrr-l];
hy3 = 2*hy2;
yc3=[yc2(end)+hy3:hy3:2*Hrr-l];
hy4 = (H-yc3(end))/4;
yc4=[yc3(end)+hy4:hy4:H];

yc = [yc1 yc2 yc3 yc4];
Ny= length(yc);

hx = l/10; % min step in x
hy = l/10; % min step in y

return
