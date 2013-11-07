
%% ------------------------------------- impose Dirichlet b.c. ------------
if wh=='g1'    % ... exact solution known, homogeneous material assumed
  [S_cur,K,A,AP,rhs_d,rhs_p] = Dirichlet_Esdo_exact(S_cur,K,A,AP,rhs_cur,...
                                   Node,Node_flagx,Node_flagy,nnode,...
                                   rho,q,r,H,L);
else           % ... exact solution NOT known
  [S_cur,K,A] = Dirichlet_Esdo_matrix(S_cur,K,A,Node_flagx,Node_flagy,nnode);
  [rhs_cur]   = Dirichlet_Esdo_rhs(rhs_cur,Node_flagx,Node_flagy,nnode,nnodeP);
end

