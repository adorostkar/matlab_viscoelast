% --------------------------------------------------------------------
% Evaluates the bilinear b.f. at point 'k' from array 'points'
% 'points' contains the coordinates of the reference triangle
% (-1,-1),(-1,3),(3,3)(3,-1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,3)   (3,3)
%     2    3
%      ----           
%     |    |         coarse_FUN(4)
%     |    |
%    1|____|4
%  (-1,-1)  (3,-1)
% --------------------------------------------------------------------


function [BF]=shape_fun_quadTH(Gauss_point,k)

      ksi  = Gauss_point(1,k);
      eta  = Gauss_point(2,k);
      
      ksim = 0.25*(3-ksi);
      ksip = 0.25*(1+ksi);
      etam = 0.25*(3-eta);
      etap = 0.25*(1+eta);
 
      BF(1,1) = ksim*etam;
      BF(1,2) = ksim*etap;
      BF(1,3) = ksip*etap;
      BF(1,4) = ksip*etam;                      

return
