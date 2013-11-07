% Visco-elasticity
% Solution of the arising linear system of equations
% C0S = C0o - 0.5*delta_t*C0

function uvp_cur = solve_lin_sysTH(S_cur,rhs_cur,A11,A11I,...
                             A12,A21,Z12,Sch,... 
                             Node,Node_flagx,Node_flagy,nnode,solver_type)

global avcgit_d avcgit_u
global kkk lll

s  = size(S_cur,1);
s1 = size(A11,1);
s2 = s - s1;

eps_gcg       = 1e-6;
max_vec       = s;
max_iter      = s;
absrel        = 'rel';
nonzero_guess = 'no ';
init_guess = zeros(s,1);

tic
if solver_type=='0', % direct solver
   uvp_cur = S_cur\rhs_cur;
   it_gcg= 0;avcgit_d=0;
else

[it_gcg,uvp_cur,resid]=gcgmr_blkf(S_cur,rhs_cur,init_guess,...
                                A11,A11I,A12,A21,Z12,Sch,...
                                max_vec,max_iter,eps_gcg,absrel,nonzero_guess); 
end
disp(['Solution time: ' num2str(toc)])
disp(['Number of GCG iterations:      ' int2str(it_gcg)])
disp(['Av. number of inner GCG iter.: ' int2str(avcgit_d)])

%
