% --------------------------------------------------------------------
% Integr_weights_quad:
% Numerical integration using Gauss quadrature formulas
% for quadrilaterals - reference element data
% --------------------------------------------------------------------
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----           
%     |    |       
%     |    |
%    1|____|4
%  (-1,-1)  (1,-1)
% --------------------------------------------------------------------

function [Gauss_point,Gauss_weight]=Integr_weights_quad

% Number of Gauss points
np = 4;
% Coordinates of Gauss points
vv=1/sqrt(3);
Gauss_point(1,1) = -vv; Gauss_point(2,1) = -vv;
Gauss_point(1,2) = -vv; Gauss_point(2,2) =  vv;
Gauss_point(1,3) =  vv; Gauss_point(2,3) =  vv;
Gauss_point(1,4) =  vv; Gauss_point(2,4) = -vv;

% Weights in Gauss points
Gauss_weight(1:np) = 1;

return
