% Compute the element stiffness matrices corresponding to the
% elasticity problem in saddle-point form
% FUND - shape functions for the displacements
% DERD - derivatives of FUND
% FUNP - shape functions for the pressure
% DERP - derivatives of FUNP
% -----> in this particular case they are the same
% --------------------------------------------------------------------
function [E_elem,B1_elem,B2_elem,M_elem,DerivD]=...
         Assm_ElSaddle_quadrTH(Gauss_point,Gauss_weight,...
	                       FUND_all,DERD_all,FUNPc,Ga,x2,y2,...
                               Coord,CoordP,wh)

      np    = 4;                    % number of points per f.e.
      nip   = 4;                    % nip = number of integration points
      dim   = 2;                    % problem dimension
      ndof  = 8;                    % ndof = np*dim
      L_elem  = zeros(ndof,ndof);   % Laplace part
      R_elem  = zeros(ndof,ndof);   % Rot part
      M_elem  = [];%zeros(ndof,ndof);   % Mass matrix for the elasticity
      E_elem  = zeros(ndof,ndof);   % Elasticity part
      Z       = zeros(np,np);       % Zero block
      B1_elem = zeros(np,np);       % Grad_x part
      B2_elem = zeros(np,np);       % Grad_y part
%      BasisP  = zeros(np,np);      % The nodal basis functions(i,np),i=1,..,np
% ------------------------------------------------------------------
% Gauss_point(2,nip)        Coord(nip,2)
% FUN(np)                   DER(2,np)          Jac(2,2)     IJac(2,2)
% Deriv(2,np)               B_cur(np,ndof)     D(np,np)
%
      for k=1:nip                         
%         FUND   = shape_fun_quad(Gauss_point,k);    % FUN(1,np)
%         DERD   = shape_der_quad(Gauss_point,k);    % DER(2,np),bilinear b.f.
%         FUNP   = shape_fun_quadTH(Gauss_point,k);

%         FUND   = FUND_all(:,k);   %shape_fun_brick(Gauss_point,k);     % (4x1)
         x3=x2'*Ga(:,k);         
         y3=y2'*Ga(:,k);         
         FUNP   = FUNPc'*[x3;y3;x3*y3;1];        
         DERD   = DERD_all(:,:,k); %shape_der_brick(Gauss_point,k);     % (2x4)
         JacD   = DERD*Coord;                     % JacD(2,2)=(2,nip)*(nip,2)
         DetD   = determinant2_m(JacD);
         IJacD  = inv(JacD);
         DerivD = IJacD*DERD;                         % (2xnp)=(2x2)*(2xnp)
	 L_elem0= DerivD'*DerivD; L_elem = [L_elem0,Z;Z,L_elem0];
	 G1_elem= DerivD(1,:)'*FUNP';              % (np,np)=(np,1)*(1,np)
	 G2_elem= DerivD(2,:)'*FUNP';              % (np,np)=(np,1)*(1,np)
	 R_elem = assm_rot_quad(DerivD);
         E_elem = E_elem + 2*DetD*Gauss_weight(k)*L_elem ... % (ndof x ndof)
	                 -   DetD*Gauss_weight(k)*R_elem;  % elasticity terms
	 B1_elem= B1_elem + DetD*Gauss_weight(k)*G1_elem;
	 B2_elem= B2_elem + DetD*Gauss_weight(k)*G2_elem;
%	 MTM    = pre_mass(FUND,ndof,dim);     
%         M_elem = M_elem  + DetD*Gauss_weight(k)*MTM; % OBS! not in sdo !!
      end
return


%
