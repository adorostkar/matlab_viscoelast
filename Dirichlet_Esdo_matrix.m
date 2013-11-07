% Imposes (homogeneous) Dirichlet b.c. according to a list of nodes Node_flag
%  (on the matrix and rhs)
% Node_flag(k,1)=1 -> Dirichlet
% Node_flag(k,1)=2 -> Neumann
% !!!!! separate-displacement ordering !!!!!
%
function [K,E,A]=...
Dirichlet_Esdo_matrix(K,E,A,Node_flagx,Node_flagy,nnode)

flx = find(Node_flagx(:,1)==1);
fly = find(Node_flagy(:,1)==1);

K(flx,:) = 0;
K(:,flx) = 0;
K(fly+nnode,:) = 0;
K(:,fly+nnode) = 0;

dd = spalloc(size(K,1),1,0);
dd(flx)=1;
dd(fly+nnode)=1;

K = K + diag(dd);

return
