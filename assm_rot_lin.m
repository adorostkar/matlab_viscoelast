% Assembly of the element stiffness matrix for the 'rot' operator
% 
%
% S(8x8) = [ S11(4x4) S12(4x4) ]
%          [ S21(4x4) S22(4x4) ]
%
% G is assumed to have the following structure

%     [df_1/dx  df_2/dx df_3/dx ]
%   G=[df_1/dy  df_2/dy df_3/dy ]
%     [    *        *       *   ]

function [S]=assm_rot_lin(G)%,V)
S=zeros(8,8);
% % pointwise
% for k=1:4,
%     for l=1:4,
%         S(k,l)     =  G(2,k)*G(2,l);  % dfk_y*dfl_y
%  	    S(k,l+4)   = -G(2,k)*G(1,l);  %-dfk_y*dfl_x
%         S(k+4,l)   = -G(1,k)*G(2,l);  %-dfl_x*dfk_y
%         S(k+4,l+4) =  G(1,k)*G(1,l);  % dfk_x*dfl_x
%     end
% end

S11 =  G(2,:)'*G(2,:);
S12 = -G(2,:)'*G(1,:); 
S21 = -G(1,:)'*G(2,:);
S22 =  G(1,:)'*G(1,:);
S = [S11 S12;S21 S22];
return