% Given: Node,Edge,Face
% asembling of the global stiffness matrix (Neumann everywhere)
% K(nnode,nnode) and the right-hand-side vector
% K is in saddle-point form; kinematic pressure intruduced
% Stable discretization (modified Taylor-Hood elements) !!!
% --> C  is the 22 matrix block 
%   ----------------------
%   |          |         |
%   |    3     |    4    |
%   |          |         |
%   |--------------------| 
%   |          |         |
%   |    2     |    1    |
%   |          |         |
%   ----------------------


function [S,K,B1,B2,C,Face_estifSo] = AssemblyV_El_quadrTH(Node,Face_Node,...
		                    Node_flagx,Node_flagy,...
                            Face_Parent,Face_estifSo,...
                            Face_eorder11,Face_eorder22,...
                            Face_flag,Face_thick,Disco,Discoef,...
                            L_char,S_char,U_char,... 
			                nnodeP,nfaceP,wh,delta_t)
global debug;

nip   = 4;         % number of integration points
dim   = 2;         %
ndof  = 8;         % ndof = nip*dim
no_children = 4;   % number of children per macroelement
I4  = eye(4);      % needed to compute the pressure b.f.
lvl = 2;           %

nface = size(Face_Node,2);   % number of finite elements
nnode = size(Node,2);   % number of nodes
K   = spalloc(2*nnode,2*nnode,ndof*2*nnode);
B1  = spalloc(nnode,nnode,ndof*nnode);
B2  = spalloc(nnode,nnodeP,ndof*nnode);
C   = spalloc(nnodeP,nnodeP,9*nnode);
nall = 2*nnode;
if(debug ~= 0)
    disp('Begin allocating memory...')
end
lengthK=nface*ndof*ndof;
KI  = zeros(lengthK,1);
KJ  = zeros(lengthK,1);
KV  = zeros(lengthK,1);
nextK = 0;

lengthB1=nface*nip*nip;     % B1 -> mju*p*phi_x
B1I = zeros(lengthB1,1);
B1J = zeros(lengthB1,1);
B1V = zeros(lengthB1,1);
nextB1 = 0;

lengthB2=nface*nip*nip;    % B2 -> mju*p*phi_y
B2I = zeros(lengthB2,1);
B2J = zeros(lengthB2,1);
B2V = zeros(lengthB2,1);
nextB2 = 0;

lengthC=nfaceP*nip*nip;     % the _22 block (mass matrix)
CI  = zeros(lengthC,1);
CJ  = zeros(lengthC,1);
CV  = zeros(lengthC,1);
nextC = 0;

if(debug ~= 0)
    disp('...end allocating memory.')
end

[Gauss_point,Gauss_weight] = Integr_weights_quad;
			  
for iface_p=1:nfaceP,     % ---> walk on the parent (pressure) elements
%    disp(['Pressure face ' int2str(iface_p)])
     local_node_listP = Face_eorder22(:,iface_p);

   CoordP(1:nip,1) = Node(1,local_node_listP)';  % CoordP(nip,2)
   CoordP(1:nip,2) = Node(2,local_node_listP)';

% - - - - - generation of the pressure element matrices:
    [C_elem,M_elem0]=Assm_quadrTH(Gauss_point,Gauss_weight,CoordP,wh);

   S_macro = zeros(18+4,18+4);
   K_macro = zeros(18,18);
   B1_macro= zeros(9,4);
   B2_macro= zeros(9,4);

   for iface_c=1:no_children
%    disp(['----------> Displacement face ' int2str(iface_c)])
       face_child       = Face_Parent(iface_c,iface_p,lvl-1);
       local_node_list  = Face_Node(:,face_child,lvl);
       local_node_list11= Face_eorder11(:,iface_p); 
       macro_node_list  = [];
       for ii=1:4,
           for jj=1:4,
               if local_node_list(ii)==local_node_listP(jj);
                  shft = ii; 
               end
           end
           for jj=1:9,
               if local_node_list11(jj)== local_node_list(ii),
                  macro_node_list=[macro_node_list,jj];
               end
           end
       end     

   local_node_list=circshift(local_node_list,[-(shft-1),0]);
   Coord(1:nip,1) = Node(1,local_node_list(1:nip))';  % Coord(nip,2)
   Coord(1:nip,2) = Node(2,local_node_list(1:nip))';

   CoordP(1:nip,1) = Node(1,local_node_listP)';  % CoordP(nip,2)
   CoordP(1:nip,2) = Node(2,local_node_listP)';

% ------------------------ determine nju, E
nju  = Discoef(1,Face_flag(face_child,1));
E    = Discoef(2,Face_flag(face_child,1))*Face_thick(face_child,1);
alfa = Discoef(3,Face_flag(face_child,1));
% ------------------------ compute the coefficients lan, mju
mju = E/(2*(1+nju))*alfa;
%lan = 2*mju*nju/(1-2*nju)*alfa;
rho = (1-2*nju)/(2*nju)*alfa;        %mju/lan*alfa;

% - - - - - generation of the element matrices:

    [El_elem,B1_elem,B2_elem,M_elem,DerivD]=...
           Assm_ElSaddle_quadrTH(Gauss_point,Gauss_weight,...
                             Coord,CoordP,wh);
% -------------- the following ordering is NOT for separate displacements
% local_node_liste = [2*local_node_list(1)-1,2*local_node_list(1),...
%                     2*local_node_list(2)-1,2*local_node_list(2),...
%                     2*local_node_list(3)-1,2*local_node_list(3),...
%                     2*local_node_list(4)-1,2*local_node_list(4)];

% -------------- the following ordering is for separate displacements     
local_node_liste = [local_node_list(1),local_node_list(2),...
                    local_node_list(3),local_node_list(4),...
                    local_node_list(1)+nnode,local_node_list(2)+nnode,...
                    local_node_list(3)+nnode,local_node_list(4)+nnode];	

macro_node_liste = [macro_node_list(1),macro_node_list(2),...
                    macro_node_list(3),macro_node_list(4),...
                    macro_node_list(1)+9,macro_node_list(2)+9,...
                    macro_node_list(3)+9,macro_node_list(4)+9];	

% Assembly of the global stiffness and mass matrices
    for ii=1:ndof,
       is = local_node_liste(ii);
       im = macro_node_liste(ii);
       for jj=1:ndof,
           js = local_node_liste(jj);
           jm = macro_node_liste(jj);
           K_macro(im,jm) = K_macro(im,jm)+ mju*El_elem(ii,jj);
%          K(is,js) = K(is,js) + mju*El_elem(ii,jj);
           nextK = nextK + 1;
           KI(nextK) = is;
           KJ(nextK) = js;
           KV(nextK) = mju*El_elem(ii,jj);
       end
    end

    for ii=1:nip,
       is = local_node_list(ii);
       im = macro_node_list(ii);
       for jj=1:nip,
          js = local_node_listP(jj);
          B1_macro(im,jj) = B1_macro(im,jj) + mju*B1_elem(ii,jj);     
%	  B1(is,js)= B1(is,js)+ mju*B1_elem(ii,jj);
          nextB1 = nextB1 + 1;
          B1I(nextB1) = is;
          B1J(nextB1) = js;
          B1V(nextB1) = mju*B1_elem(ii,jj);

          B2_macro(im,jj) =B2_macro(im,jj) + mju*B2_elem(ii,jj);           
%	  B2(is,js)= B2(is,js)+ mju*B2_elem(ii,jj);
          nextB2 = nextB2 + 1;
          B2I(nextB2) = is;
          B2J(nextB2) = js;
          B2V(nextB2) = mju*B2_elem(ii,jj);
 
       end
    end

   end   % end children's face loop
   
   for ii=1:nip,
       is = local_node_listP(ii);
       for jj=1:nip,
          js = local_node_listP(jj);

%         C(is,js) = C(is,js) - rho*mju*M_elem0(ii,jj);
         nextC = nextC + 1;
         CI(nextC) = is;
         CJ(nextC) = js;
         CV(nextC) = - rho*mju*M_elem0(ii,jj);     
       end
    end
    C_macro = - rho*mju*M_elem0;
    B_macro = [B1_macro;B2_macro];
    S_macro = [K_macro B_macro; B_macro' C_macro];
    Face_estifSo(1:22,1:22,iface_p) = S_macro;
    
end             % end pressure iface_p-loop


K  = sparse(KI, KJ, KV);
B1 = sparse(B1I,B1J,B1V);
B2 = sparse(B2I,B2J,B2V);
C  = sparse(CI, CJ, CV); 

B = [B1;B2];
S = [K B; B' C];

return


    
