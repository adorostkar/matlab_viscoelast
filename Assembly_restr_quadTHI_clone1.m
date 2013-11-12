% This program assembles element Schur A11I and Z12 element blocks

% ndofe = numder of degree of freedom per finite element
% ndofm = number of degrees of freedom per element
% ndofc = number of coarse nodes per element
% nface = number of faces 
%
% The block A11 corresponds to the matrix 'S - theta*delta_t_cur*So'
% and all has to be recomputed if theta or delta_t are changed !!!
%

function [A11I,Z12,Z21,S,...
         Face_eSchur,Face_e11inv]=...
         Assembly_restr_quadTHI(A11as,Node,...
                Face_estiffS,Face_estiffSo,Face_eorder11,Face_eorder22,...
                Face_eSchur,Face_e11inv,...
                Node_Face_no,nface,ndoff,ndofc,delta_t,theta);

global h

nnode = size(Node,2);
ndofm = ndoff+ndofc;

lengthS = nface*ndofc*ndofc;
SI  = zeros(lengthS,1);
SJ  = zeros(lengthS,1);
SV  = zeros(lengthS,1);
nextS   = 0;

%lenA11I = nface*ndoff*ndoff;
%A11II = zeros(lenA11I,1);
%A11IJ = zeros(lenA11I,1);
%A11IV = zeros(lenA11I,1);
%nextA11I = 0;


%lenZ12 = nface*ndoff*ndofc;
%Z12I = zeros(lenZ12,1);
%Z12J = zeros(lenZ12,1);
%Z12V = zeros(lenZ12,1);
%nextZ12 = 0;

%lenZ21 = nface*ndoff*ndofc;
%Z21I = zeros(lenZ21,1);
%Z21J = zeros(lenZ21,1);
%Z21V = zeros(lenZ21,1);
%nextZ21 = 0;

DI = eye(ndoff);
DP = h^2*DI;         % diagonal perturbation to A11
DDS= zeros(ndoff);
DD = DI;
for iface=1:nface,   % ---> walk on faces

% ---> Load the element matrix    
  S_elem = Face_estiffS(:,:,iface);
  So_elem= Face_estiffSo(:,:,iface);

  local_node_list  = Face_eorder11(:,iface);
  local_node_listp = Face_eorder22(:,iface);
% -------------- the following ordering is NOT for separate displacements
% local_node_liste = [2*local_node_list(1)-1,2*local_node_list(1),...
%                     2*local_node_list(2)-1,2*local_node_list(2),...
%                     2*local_node_list(3)-1,2*local_node_list(3),...
%                     2*local_node_list(4)-1,2*local_node_list(4)];

% -------------- the following ordering is for separate displacements     
 local_node_liste = [local_node_list;local_node_list+nnode];	


  A_elem = S_elem - theta*delta_t*So_elem;

% Assembly the Schur complement matrix

%  for r=1:ndofc,
%      DD(r,r)  = DI(r,r)/Node_Face_no(local_node_liste(r));  
%      DD(r+ndofc,r+ndofc) = DD(r,r);
%  end

% --------- Element matrix blocks computations  
 A11 = A_elem(1:ndoff,1:ndoff);
 A12 = A_elem(1:ndoff,ndoff+1:ndofm);
 A21 = A_elem(ndoff+1:ndofm,1:ndoff);
 A22 = A_elem(ndoff+1:ndofm,ndoff+1:ndofm);

% ---> Restrict A11 to the degrees of freedom per element  
%      from the assembled A11 block (A11as)
 A11R =  A11as(local_node_liste,local_node_liste);
% !! NO scaling for the Schur complement  

 invA11    = inv(A11R);% ---------->  OOBBSS ! Two inverses !!!
 invA11S   = inv(A11+DP); %              How to avoid them ????  
 Sch_elem  = A22 - A21*invA11S*A12;

 Face_eSchur(1:ndofc,1:ndofc,iface) = Sch_elem;
%  invA11 = invA11*DD;
%  Z12_elem = invA11S*DD*A12; 
%  Z21_elem = A21*invA11S*DD;
%  Face_e11inv(1:ndoff,1:ndoff,iface) = invA11;
%  Face_e11x12(1:ndoff,1:ndofc,iface) = Z12_elem;
%  Face_e21x11(1:ndofc,1:ndoff,iface) = Z21_elem;

   for i=1:ndofc,
       is = local_node_listp(i);
      for j=1:ndofc,
         js = local_node_listp(j);
%         S(is,js) = S(is,js) + Sch_elem(i,j);          
         nextS = nextS + 1;
         SI(nextS) = is;
         SJ(nextS) = js;
         SV(nextS) = Sch_elem(i,j); 
      end
   end

%    for i=1:ndoff,
%        is = local_node_liste(i);
%       for j=1:ndoff,
%          js = local_node_liste(j);
%%         A11I(is,js) = A11I(is,js) + invA11(i,j);          
%          nextA11I = nextA11I + 1;
%          A11II(nextA11I) = is;
%          A11IJ(nextA11I) = js;
%          A11IV(nextA11I) = invA11(i,j); 
%       end
%    end
	  
%    for i=1:ndoff,
%        is = local_node_liste(i);
%       for j=1:ndofc,
%          js = local_node_listp(j);
%%         Z12(is,js) = Z12(is,js) + Z12_elem(i,j); 
%          nextZ12 = nextZ12 + 1;
%          Z12I(nextZ12) = is;
%          Z12J(nextZ12) = js;
%          Z12V(nextZ12) = Z12_elem(i,j);
%	          
%%         Z21(is,js) = Z21(is,js) + Z21_elem(j,i); 
%          nextZ21 = nextZ21+ 1;
%          Z21I(nextZ12) = js;
%          Z21J(nextZ12) = is;
%          Z21V(nextZ12) = Z21_elem(j,i);	          
%       end
%    end


end    % ---> end walk on faces

S    = sparse(SI, SJ, SV);          %
A11I = [];%sparse(A11II, A11IJ, A11IV); %
Z12  = [];%sparse(Z12I, Z12J, Z12V);    %
Z21  = [];%sparse(Z21I, Z21J, Z21V);    %

return