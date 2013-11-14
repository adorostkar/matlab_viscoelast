% Numerical integration of a function f
% along an line interval, given by two points: (x1,y1);(x2,y2)
% The integration formula involves the midpoint of the interval
% The load can be time-dependent
%

function [F,len] = lin_num_int2L_var(x1,y1,x2,y2,...
                      wh,l_ice,h_ice,rho_ice,grav,...
                      time_length,T_BEG,T_LGM,T_EOG)

xdif = abs(x1-x2);
len  = sqrt((x2-x1)^2+(y2-y1)^2);  % the length of the edge
cosa = xdif/len;
% surface_force_var should compute the force/pressure due to ice load
% in the moddle of the edge, incorporating the edge does not have to be
% horizontal
F  = surface_force_var(x1,y1,x2,y2,cosa,...
                        wh,l_ice,h_ice,rho_ice,grav,...
                        time_length,T_BEG,T_LGM,T_EOG);
F  = 0.5*F*len;
return 