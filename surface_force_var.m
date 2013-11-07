% Gives the value of the surface force
% wh is a problem identifier
% ---------------------------------------------------------
% wh='gs': vertical load originally on 0<=x<=l_ice, y=0
%          horizontal load nowhere
% ---------------------------------------------------------
% wh='g1': vertical load   nowhere
%          horizontal load nowhere
% ---------------------------------------------------------
% wh='g2': vertical   load (changing linearly)
%          horizontal load nowhere
% ---------------------------------------------------------
% OBS! When compuing the surface force we multiply by 0.5
%      (which is he value of the basis function in the middle
%       of the interval)
% ---------------------------------------------------------
% OBS! The function is called ONLY for the edges under the ice,
%      so, the chech where the ice is no longer needed

function F = surface_force_var(x1,y1,x2,y2,cosa,...
                     wh,l_ice,h_ice,rho_ice,grav,...
                     time_length,T_BEG,T_LGM,T_EOG)

F(1) = 0; % horizontal force
F(2) = 0; % vertical   force

if   (wh=='g1'), % no surface forces,  quadrilaterals
    return
else
% determine the current height of the ice at the point (x,y)
     x = (x1+x2)*0.5;
     y = (y1+y2)*0.5;
     h_ice_cur = height_ice(wh,x,y,h_ice,time_length,T_BEG,T_LGM,T_EOG);
     F(2) = -grav*rho_ice*h_ice_cur*0.5*cosa;
end

return
