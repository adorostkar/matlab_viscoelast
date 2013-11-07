% --------------------------------------------------------------------
%
% Evaluates the first partial derivatives of the bilinear b.f.
% at point 'k' from array 'Gauss_point'
% 'Gauss_point' contains the coordinates of the
%  reference triangle (0,0),(1,0),(0,1)
%
% Presummed ordering of the nodepoints in the reference f.e.:
%  (-1,1) (1,1)
%     2    3
%      ----           
%     |    |       DER = [ dfi_i/dksi; dfi_i/deta]
%     |    |       DER(2,4)
%    1|____|4
%  (-1,-1)  (1,-1)

%      BF(1,1) = 4*ksim*etam;
%      BF(1,2) = 4*ksim*etap;
%      BF(1,3) = 4*ksip*etap;
%      BF(1,4) = 4*ksip*etam;
% --------------------------------------------------------------------


function [DER]=shape_der_quad(Gauss_point,k) %

      ksi  = Gauss_point(1,k);
      eta  = Gauss_point(2,k);
      ksim = 0.25*(1-ksi);
      ksip = 0.25*(1+ksi);
      etam = 0.25*(1-eta);
      etap = 0.25*(1+eta);

      DER(1,1) = -etam;
      DER(2,1) = -ksim;

      DER(1,2) = -etap;
      DER(2,2) =  ksim;

      DER(1,3) =  etap;
      DER(2,3) =  ksip;

      DER(1,4) =  etam;
      DER(2,4) = -ksip;

return
