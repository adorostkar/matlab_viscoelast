% Computes the current height of the ice cap as a function of time
% 

function h = height_ice(wh,x,y,h_ice,time_length,T_BEG,T_LGM,T_EOG)

if (wh=='gs')|(wh=='g0'),
   h = h_ice; % boxcar with fixed height
% currentlyy all other ice shapes are excluded/not tested
elseif wh=='g2',
     T0 = T_LGM-T_BEG;
     T1 = T_EOG-T_LGM;
     if time_length<T_BEG,
       h = 0;
     elseif time_length<T_LGM,
       h  = h_ice/T0*(time_length- T_BEG);
     elseif time_length<T_EOG,
       h  = h_ice/T1*(time_length- T_EOG);
     else
       h = 0;
     end   
elseif wh=='g1',
     h = 0;
else
   disp('visco_elast_rhs: Wrong problem identifier.'),return     
end

return
