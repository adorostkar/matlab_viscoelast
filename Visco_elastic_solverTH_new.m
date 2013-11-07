%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The viscoelastic solver
% ----- Initialize:
% I1. Assembly the stiffness matrix,  S = A_elast + A_advec,
%     elasticity material coefficients at time tau=t
% I2. Assembly the elasticity stiffness matrix A_visc(0)
%     viscoelast coeff at time tau=t
% I3. Set k=1 add choose time t(k) and timestep dt(k)=t(k)-t(k-1)
% I4. Compute the initial rhs0
% I5. Choose initial guess for displacements and pressure, U(0)
% I6. Assembly the elasticity stiffness matrix A_visc(k)
%     viscoelast. coeff at time tau=dt(k)
% I7. W(k) = dt(k)/2*A_visc(k)*U(0)
%     rhs = rhs0 + W(0)
% I8. Solve (S - dk/2*A_visc(k))*U(k) = rhs;
% I9. Compute the integral [0->t(k)]:
%     I(K) = dk/2*A_visc(0)*U(k) + W(k)
% ----- Loop over time
% L1. Set k=k+1, determine timestep dt(k) = t(k)-t(k-1)
% L2. Assembly the elasticity stiffness matrix A_visc(k)
%     viscoelast coeff at time tau=tau+dk
% L3: Compute the rhs corresponding to a new load - rhs(k)
% L4. W(k) = dt(k)/2*A_visc(k)*U(k-1)
%     rhs  = rhs(k) + W(k) + exp(-dk/Maxwell_time)*I(k)
% L5. Solve (S - dt(k)/2*A_visc0)*U(k) = rhs;
% L6. Compute the integral term:
%     I(k)=I(k-1) + 1/2*dk(k-1)*A_visc(k)*U(k) + W(k)
% L7. Goto L1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I1. Assembly the global stiffness matrix
%     elasticity material coefficients at time tau=t
% K - elasticity part of S11
% A - advection part of S11  (S11=K+A)
%     [K+A     [B1 B2]^T ]
% S = [[B1 B2]    C      ]
% K and A - in separate displacement ordering
% C0 - pressure mass matix

%ilu_setup.droptol = actionILU(kkk);

T_BEG = 0; % the time point from which the ise load is imposed

disp('WITH PRE-STRESS ADVECTION.')  % with advection
vec_coeff0 = [0,-rho_earth*grav,0,-rho_earth*grav];
%vec_coeff0 = [0,-0.5,0,-0.5],pause
%disp('NO PRE-STRESS ADVECTION.')  % no advection
%vec_coeff0 = [0,0,0,0];
theta = 0.5;    % in this case we use Trapetz method

vec_coeff  = L_char/S_char*vec_coeff0;

Face_estiffS =zeros(ndof+np,ndof+np,nface);

[S,K,A,AP,B1,B2,...
 C,Face_estiffS,...
 rhs_db,rhs_pb] = AssemblyV_ElAd_quadrTH(Node,Face_Node,...
                           Node_flagx,Node_flagy,...
                           Face_Parent,Face_estiffS,...
                           Face_eorder11,Face_eorder22,...
                           Face_flag,Face_thick,Disco,Discoef,...
                           L_char,S_char,U_char,...
                           nnodeP,nfaceP,vec_coeff,wh); %pure elastic


% I2. Assembly the elasticity stiffness matrix A_visc0
%     viscoelast coeff at time tau=t
Face_estiffSo =zeros(ndof+np,ndof+np,nface);
[So,Ko,B1o,B2o,...
 Co,Face_estiffSo] = AssemblyV_El_quadrTH(Node,Face_Node,...
                              Node_flagx,Node_flagy,...
                              Face_Parent,Face_estiffSo,...
                              Face_eorder11,Face_eorder22,...
                              Face_flag,Face_thick,Disco,Discoef,...
                              L_char,S_char,U_char,...
                              nnodeP,nfaceP,wh,0);  % delta_t = 0

% I3. Choose initial guess for displacements and pressure (zero)
Maxwell_time_inv = Discoef(3,1);%
% Time scaling: spY     = 365*24*60*60; - seconds per year
%               Years   = 200000;
%               T0      = Years*spY
%               T_char  = N_char/S_char;
%               T       = T0/T_char;
%               delta_tY= delta_t*T_char/spY (the timestep in years)

%delta_t          = 0.01;% 0.01 = 0.0113 years
delta_t          = 0.87;  %  0.87 = 1 year
% delta_t          = 69.60; % 69.60 = 80 years
delta_t_cur      = delta_t;
time_length      = 0;
time_length      = time_length + delta_t_cur;

% I4. Compute the initial rhs0
[rhs_ds,rhs_ps] = Assembly_rhsTH(Node,Edge,...
                          wh,l_ice,h_ice,rho_ice,grav,Bndry_Ice,...
                          L_char,S_char,U_char,nnodeP,...
                          time_length,T_BEG,T_LGM,T_EOG);

rhs_d   = rhs_db+rhs_ds;
rhs_p   = rhs_pb+rhs_ps;
rhs_cur = [rhs_d;rhs_p];

% I5. Choose initial guess for displacements and pressure (zero)
uvp_prev        = zeros(2*nnode+nnodeP,1);
UVPx(1:nnode, 1) = zeros(nnode,1);
UVPy(1:nnode, 1) = zeros(nnode,1);
UVPp(1:nnodeP,1) = zeros(nnodeP,1);

% I6. Assembly the elasticity stiffness matrix A_visc(k)
%     viscoelast coeff at time tau=dk
Face_estiffSo =zeros(ndof+np,ndof+np,nface);
[Sc,Kc,B1c,B2c,...
 Cc,Face_estiffSo] = AssemblyV_El_quadrTH(Node,Face_Node,...
                              Node_flagx,Node_flagy,...
                              Face_Parent,Face_estiffSo,...
                              Face_eorder11,Face_eorder22,...
                              Face_flag,Face_thick,Disco,Discoef,...
                              L_char,S_char,U_char,...
                              nnodeP,nfaceP,wh,delta_t);

% I7. W(k) = dt(k)/2*A_visc(k)*U(0); rhs = rhs0 + W(0)
W = theta*delta_t_cur*Sc*uvp_prev;   % W = 0
rhs_cur = rhs_cur + W;               % kept for consistency

% I8. Solve (S - dk/2*A_visc(k))*U(k) = rhs;

S_cur = S - theta*delta_t_cur*So;
%S_cur = S;
Dirichlet_visco            %  impose b.c. both to S_cur and rhs_cur

%  compute U_1 = U(t1)
uvp_cur = S_cur\rhs_cur;   %  compute U_1 = U(t1)

% --------------------------------    Plot the current solution
UVPx(1:nnode, 2) = uvp_cur(1:nnode,1);
UVPy(1:nnode, 2) = uvp_cur(nnode+1:2*nnode,1);
UVPp(1:nnodeP,2) = uvp_cur(2*nnode+1:2*nnode+nnodeP,1);

% I9. Compute the integral [0->t(k)]:
%     I(K) = dk/2*A_visc(0)*U(k) + W(k)
I = delta_t_cur*theta*So*uvp_cur + W;

% ------------------------------- Loop over time
k = 1;
test_regime=0;
norm_dif = norm(uvp_cur-uvp_prev);
norm_cur = norm(uvp_cur);
norm_prev = norm(uvp_prev);
% while (norm(uvp_cur-uvp_prev)>1e-6)&(time_length<=Tmax)
% while (time_length<=Tmax && k < 4)
while (norm_dif>1e-6)&(time_length<=Tmax)
    k = k + 1;
    
    fprintf('Updating Node placements\n');
    for i=1:nnode
        Node(1,i) = UVPx(i,k) + Node(1,i);
        Node(2,i) = UVPy(i,k) + Node(2,i);
    end
   if(verbose ~= 0)
       figure(2),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,0,16)
   end
 
    %     figure
    %     plot(Node(1,Surface_Nodes), UVPx(Surface_Nodes,end),'.')
    
    disp(['Proceed with step ' int2str(k)]),
    disp(['norm(uvp_cur-uvp_prev) ' num2str(norm_dif) ...
        ' ,norm(uvp_cur) ' num2str(norm_cur) ...
        ' ,norm(uvp_prev) ' num2str(norm_prev)])
    
    [S,K,A,AP,B1,B2,C,...
     Face_estiffS,...
     rhs_db,rhs_pb] = AssemblyV_ElAd_quadrTH(Node,Face_Node,...
                               Node_flagx,Node_flagy,...
                               Face_Parent,Face_estiffS,...
                               Face_eorder11,Face_eorder22,...
                               Face_flag,Face_thick,Disco,Discoef,...
                               L_char,S_char,U_char,...
                               nnodeP,nfaceP,vec_coeff,wh); %pure elastic
    
    [So,Ko,B1o,B2o,...
     Co,Face_estiffSo] = AssemblyV_El_quadrTH(Node,Face_Node,...
                                  Node_flagx,Node_flagy,...
                                  Face_Parent,Face_estiffSo,...
                                  Face_eorder11,Face_eorder22,...
                                  Face_flag,Face_thick,Disco,Discoef,...
                                  L_char,S_char,U_char,...
                                  nnodeP,nfaceP,wh,0);  % delta_t = 0
    
    delta_t_prev = delta_t_cur;
    uvp_prev     = uvp_cur;
    
    % L1. Set k=k+1, determine timestep dt(k) = t(k)-t(k-1)
    delta_t_cur = delta_t_prev;  % take the same time step
    time_length = time_length + delta_t_cur;
    
    % L2. Assembly the elasticity stiffness matrix A_visc(k)
    %     viscoelast coeff at time tau=tau+dk
    if (delta_t_prev~=delta_t_cur)
        [Sc,Kc,B1c,B2c,...
         Cc,C0c] = AssemblyV_El_quadrTH(Node,Face_Node,...
                               Node_flagx,Node_flagy,...
                               Face_Parent,Face_estiffS,...
                               Face_eorder11,Face_eorder22,...
                               Face_flag,Face_thick,Disco,Discoef,...
                               L_char,S_char,U_char,...
                               nnodeP,nfaceP,vec_coeff,wh);
    end
    
    % L3: Compute the rhs corresponding to a new load - rhs(k)
    
    [rhs_ds,rhs_ps] = Assembly_rhsTH(Node,Edge,...
                              wh,l_ice,h_ice,rho_ice,grav,Bndry_Ice,...
                              L_char,S_char,U_char,nnodeP,...
                              time_length,T_BEG,T_LGM,T_EOG);

    rhs_d   = rhs_db+rhs_ds;
    rhs_p   = rhs_pb+rhs_ps;
    rhs_cur = [rhs_d;rhs_p];
    
    % L4. W(k) = dt(k)/2*A_visc(k)*U(k-1)
    %     rhs  = rhs(k) + W(k) + exp(-dk/Maxwell_time)*I(k)
    W            = theta*delta_t_cur*Sc*uvp_prev;
    rel_time_cur = -Maxwell_time_inv*delta_t_cur;
    I            = exp(rel_time_cur)*I;
    rhs_cur      = rhs_cur + W + I;
    rhs_cur      = Dirichlet_Esdo_rhs(rhs_cur, ...
                            Node_flagx,Node_flagy,nnode,nnodeP);
    
    % L5. Solve (S - dt(k)/2*A_visc0)*U(k) = rhs;
    if  (delta_t_prev~=delta_t_cur)
        S_cur       = S - theta*delta_t_cur*So;
        [S_cur,K,A] = Dirichlet_visco_matrix(S_cur,K,A,Node,...
                               Node_flagx,Node_flagy,nnode);
    end
    
    uvp_cur = S_cur\rhs_cur;
    
    
    UVPx(1:nnode, k+1)=uvp_cur(1:nnode,1);
    UVPy(1:nnode, k+1)=uvp_cur(nnode+1:2*nnode,1);
    UVPp(1:nnodeP,k+1)=uvp_cur(2*nnode+1:2*nnode+nnodeP,1);
    
    % L6. Compute the integral term:
    %     I(k)=I(k-1) + 1/2*dk(k-1)*A_visc(k)*U(k) + W(k)
    I = I + W + theta*delta_t_cur*So*uvp_cur;
    
    % Update norms
    norm_dif = norm(uvp_cur-uvp_prev);
    norm_cur = norm(uvp_cur);
    norm_prev = norm(uvp_prev);
    
    test_regime=1; %wait
end   % L7. Goto L1

disp(['norm(uvp_cur-uvp_prev) ' num2str(norm(uvp_cur-uvp_prev)) ...
    ' ,norm(uvp_cur) ' num2str(norm(uvp_cur))])

return
