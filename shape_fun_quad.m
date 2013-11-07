% --------------------------------------------------------------------
% Evaluates the bilinear b.f. at point 'k' from array 'points'
% 'points' contains the coordinates of the reference square
% (-1,-1),(-1,1),(1,1)(1,-1),
%
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----           
%     |    |         FUN(4)
%     |    |
%    1|____|4
%  (-1,-1)  (1,-1)
% --------------------------------------------------------------------


function [BF]=shape_fun_quad(Gauss_point,k)

      ksi  = Gauss_point(1,k);
      eta  = Gauss_point(2,k);
      
      ksim = 0.25*(1-ksi);
      ksip = 0.25*(1+ksi);
      etam = 0.25*(1-eta);
      etap = 0.25*(1+eta);

      BF(1,1) = 4*ksim*etam;
      BF(1,2) = 4*ksim*etap;
      BF(1,3) = 4*ksip*etap;
      BF(1,4) = 4*ksip*etam;

return
