% Initial quadrilateral mesh
% b.c. such that the solution is displ_x=0; displ_y - linear
%
% Face_flag(iface)=disco_region
%---->  Regular meshsize in each direction
% Data: Nx - number of intervals in x-direction
%       Ny - number of intervals in y-direction
%       [hx] - variable stepsize in x-direction
%       [hy] - variable stepsize in y-direction

function [xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_my(L,H,Nx,Ny)
% ------- subdomain 1:
% OBS !!! H is given with its true sign

hx = L/(Nx-1); 
xc =[0:hx:L];  

hy = H/(Ny-1);
yc  =[0:hy:H];

return
