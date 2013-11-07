% Application of a block-factorized precinditioner
% with approximated inv(A11), inv(A11)*A12 and approx Schur complement S
% [ A11   0 ]
% [         ]
% [ A21   S ]
%
%  1: x1 = A11I*b1;        
%     x2 = S\(b2-A21*x1)    
%     (exact solve with S)
%

% ---------------------------------------------------------------
function [x,itrd,itru]=tril_prec(b,A11,A11I,A21,S,M1,M2)
global actionILU actionFEM
global L11d U11d L11u U11u% lu(A11)/luinc(A11)
global kkk lll
global eps_inner
global p_amd

itrd = 0; itru = 0;
[s2,s1]=size(A21);
s  = s1+s2;
x10=[];
b1 = b(1:s1,1);
b2 = b(s1+1:s,1);
if kkk==1,
   if lll==1,   % exact solve with the _11 block
      bb1 = b1(p_amd);
      xx1 = A11\bb1;
      x1(p_amd,1) = xx1;
   else         % block-diag gcg for the _11 block
      [it_blkdiag,x1,resid_bd]=gcgmr_blockdiag(A11,b1,x10,M1,M2,...
                           s1,s1,eps_inner,'rel','no ');
      itrd = itrd + it_blkdiag;       
   end
elseif kkk<5,
   z1 = L11d\b1;
   x1 = U11d\z1;
elseif kkk==5  % do not solve with A11, but multiply by the approx. inv
   x1 = A11I*b1;
else % kkk>5,
   maxitd = actionFEM(kkk-4);
%   [itrd,y1,resid]=gcgmr_mult(A11,b1,b1,A11I,s1,maxitd,1e-3,'rel','no ');
   [itrd,y1,resid]=gcgmr_mult(A11,b1,b1,A11I,s1,maxitd,1e-5,'rel','no ');
   disp(['inner iter.(down):  ' int2str(itrd)])
end
bb = b2-A21*x1;
x2 = S\bb;      % only direct solve with the Schur complement

x=[x1;x2];

return
% ---------------------------------------------------------------

	 

     
