% Compute the element mass matrices 
% FUND - shape functions for the displacements
% DERD - derivatives of FUND
% FUNP - shape functions for the pressure
% DERP - derivatives of FUNP
% -----> in this particular case they are the same
% coeff - ???? add description
% --------------------------------------------------------------------
function [C_elem,M_elem]=...
         Assm_quadrTH(Gauss_point,Gauss_weight,Coord,wh)

      np    = 4;                    % number of points per f.e.
      nip   = 4;                    % nip = number of integration points
      dim   = 2;                    % problem dimension
      ndof  = 8;                    % ndof = np*dim
      C_elem  = zeros(np,np);       % Laplace
      M_elem  = zeros(np,np);       % Mass matrix for the Laplacian
% ------------------------------------------------------------------
% Gauss_point(2,nip)        Coord(nip,2)
% FUN(np)                   DER(2,np)          Jac(2,2)     IJac(2,2)
% Deriv(2,np)               B_cur(np,ndof)     D(np,np)
%
      for k=1:nip                         
         FUN    = shape_fun_quad(Gauss_point,k);    % FUN(1,np)
         DER    = shape_der_quad(Gauss_point,k);    % DER(2,np),bilinear b.f.
	 Jac    = DER *Coord;                     % JacD(2,2)=(2,nip)*(nip,2)
         Det    = determinant2_m(Jac);
         IJac   = inv(Jac);
         Deriv  = IJac*DER;                         % (2xnp)=(2x2)*(2xnp)
	 C      = Deriv'*Deriv;
	 M      = FUN'*FUN;
         C_elem = C_elem + Det*Gauss_weight(k)*C;   % pressure Laplace
         M_elem = M_elem + Det*Gauss_weight(k)*M;   % pressure mass matrix
      end
return


%
