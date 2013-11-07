% Solve the inner  system H*(tau)=b,
% where tau are coeff. for: x(new)=x(old)+sum[tau(j)*d(j)]
% Parameters
% r=current residual; 

function [tau,H,flag]=solvesubsys(r,Ad,H,adlast,j0,it,max_vec,flag);

dimH = it-j0+1;   % the dimension of the small system H
if it>max_vec,
   for i=1:max_vec-1,
       for j=1:max_vec-1,
           H(i,j) = H(i+1,j+1);
       end
   end
end
%
for j=j0:it
    ThisPos        = mod(j-1,max_vec) + 1;
    H(j-j0+1,dimH) = Ad(:,ThisPos)'*adlast; H(dimH,j-j0+1) = H(j-j0+1,dimH);
    b(j-j0+1,1)    = 0;
end
b(dimH,1) = -r'*adlast;
tau = H\b;

condH=rcond(H);
if condH < 1e-12, flag = flag-1; end   % H is illconditioned, stagnaton 
%sadlast,wait
return



