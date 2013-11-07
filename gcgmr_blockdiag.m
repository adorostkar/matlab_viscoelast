% The preconditioned Generalized Conjugate Gradient - minimal residual
% implementation ;  mltipliative preconditioner
% Parameters: y0-zero initial guess; rhs-right_hand_side
% 

function [it,x,resid]=gcgmr_blockdiag(A,rhs,x,M1,M2,max_vec,max_iter,... 
                           eps_gcgmr,absrel,nonzero_guess)                          

resid= [];Hsub=[];
flagH= 0;
it   = 0;
s1 = size(M1,1);
s2 = size(M2,1);
if nnz(rhs)==0, it = 0; x=zeros(size(rhs)); resid=0; return,end
if nonzero_guess=='yes', 
   r = A*x-rhs;
else
   x = zeros(length(rhs),1);
   r = -rhs;
end 
r1=r(1:s1,1);
r2=r(s1+1:s1+s2,1);  
h1=M1\r1; h2=M2\r2; h=[h1;h2];
h=-h;
InitRes=sqrt(r'*r);resid=[resid;InitRes];
%disp(' ')
%disp(['_________________________gcgmr_rnorm0:',num2str(InitRes)])
if absrel=='rel',
   ceps = eps_gcgmr*InitRes;
else   
   ceps = eps_gcgmr;
end 
 
rnorm=1e13;
    while (rnorm>ceps)&(it<max_iter),
        it = it + 1;
        ThisPos = mod(it-1, max_vec) + 1;
        d(:,ThisPos)=h;
        Ad(:,ThisPos)=A*h;
        if it<=max_vec,
           j0 = 1;
        else
           j0 = it - max_vec + 1;
        end
       [tau,Hsub,flagH]=solvesubsys(r,Ad,Hsub,Ad(:,ThisPos),j0,it,max_vec,flagH);
        for j=j0:it,
            ThisPos = mod(j-1,max_vec) + 1;
            x = x + tau(j-j0+1)*d(:,ThisPos);
            r = r + tau(j-j0+1)*Ad(:,ThisPos);
        end
%        figure(1),plot(x),figure(2),plot(r),figure(3),plot(A*x-rhs),wait
        rnorm=sqrt(r'*r);
%        disp(['_________________________gcgmr_rnorm :',num2str(rnorm)])
        r1=r(1:s1,1);
        r2=r(s1+1:s1+s2,1);  
        h1=M1\r1; h2=M2\r2; h=[h1;h2];
        h = -h;
        resid=[resid;rnorm];
    end
return
	 

     
