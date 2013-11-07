% The preconditioned Generalized Conjugate Gradient - minimal residual
% implementation  
% Parameters: y0-zero initial guess; rhs-right_hand_side
% 

function [it,x,resid]=gcgmr_blkf(A,rhs,x,...
                           A11,A11I,A12,A21,A12t,S,max_vec,max_iter,... 
                           eps_gcgmr,absrel,nonzero_guess)
global avcgit_d avcgit_u
global fid
global p_amd
global kkk lll

s1 = size(A11,1); sh = s1/2;
M1 = A11(1:sh,1:sh);
M2 = A11(sh+1:s1,sh+1:s1);

if lll==1,
   p_amd = symamd(A11);
   A11 = A11(p_amd,p_amd);
end   
if max_iter ==0, max_iter = size(A,1); end
resid=[];Hsub=[]; 
avcgit_d = 0; avcgit_u = 0; flagH = 0;
it=0;
if nonzero_guess=='yes', 
   r = A*x-rhs;
else
   x = zeros(length(rhs),1);
   r = -rhs;
end  
 
%[h,inner_itd,inner_itu]=blkf_prec(r,A11,A11I,A12,A21,A12t,S);
[h,inner_itd,inner_itu]=tril_prec(r,A11,A11I,A21,S,M1,M2);
avcgit_d=avcgit_d+inner_itd;
avcgit_u=avcgit_u+inner_itu;
h=-h;
InitRes=sqrt(r'*r);resid=[resid;InitRes];
disp(' ')
disp(['_________________________gcgmr_rnorm0:',num2str(InitRes)])
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
	if flagH==-10, 
           max_iter= 0; 
           disp('STAGNATION!!')
           fprintf(fid,' ---STAGNATION at iteration %i \n',it);
           end
        for j=j0:it,
            ThisPos = mod(j-1,max_vec) + 1;
            x = x + tau(j-j0+1)*d(:,ThisPos);
            r = r + tau(j-j0+1)*Ad(:,ThisPos);
        end
%        figure(1),plot(x),figure(2),plot(r),figure(3),plot(A*x-rhs),wait
        rnorm=sqrt(r'*r);
        disp(['_________________________gcgmr_rnorm :',num2str(rnorm)])
%        [h,inner_itd,inner_itu]=blkf_prec(r,A11,A11I,A12,A21,A12t,S);
        [h,inner_itd,inner_itu]=tril_prec(r,A11,A11I,A21,S,M1,M2);
        avcgit_d=avcgit_d+inner_itd;
        avcgit_u=avcgit_u+inner_itu;
        h = -h;
        resid=[resid;rnorm];
    end
    avcgit_d=ceil(avcgit_d/(it+1));
    avcgit_u=ceil(avcgit_u/(it+1));
return
	 

     
