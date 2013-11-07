%               Modified Taylor-Hood stable discretization 
%
%%%%%%%%%%%%%%%%%% 2D elasticity with advection terms %%%%%%%%%%%%%%%%%%%%%%%%%%
% Works on any given coarse mesh, described as {Node, Edge, Face}
% 1. The coarse mesh is refined 'levels' times
% 2. The mesh obtained in 1. is the pressure mesh. It is then refined
%    once more to get the mesh for the displacements. At this level
%    parent-child relation is kept to utilize the assembly of the 
%    stiffness matrices 
% 3. The viscoelastic solver is run It includes the assembly of the 
%    corresponding STIFFNESS and MASS matrices 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Block-factorized preconditioner, two-by-two block structure introduced
% Parent-child relation in the mesh refinement introduced
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% ------------------------------- global definitions
global average_inner how_many outer_solv inner_solv inner_solver_type
global restart restart_tol
global avcgit_d avcgit_u
global sigma
global h
global solver_type
global actionILU actionFEM
global L11d U11d L11u U11u % lu(A11)/luinc(A11)
global kkk lll
global eps_inner
global theta

verbose = 1;
np    = 4;  % number of points in the finite element
nip   = 4;  % number of integration points
dim   = 2;  % number of degree of freedom (elasticity part)
ndof  = 8;  % ndof = nip*dim

actionILU =[999,1e-1,1e-2,1e-3];
actionILUt=['exact solve(A11)','cholinc(0.0100) ','cholinc(0.0010) ','cholinc(0.0001) '];
actionFEM =[0,1,2,3,4,5,1e6,999];
actionFEMt=['iter=0 ','iter=1 ','iter=2 ','iter=3 ','iter=4 ',...
            'iter=5 ','iter=50','Z12    '];
lll = 2;
kkk = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Homogeneous half-space %%%%%%%%%%%%
% wh = 'gs' 
% Nx = 10;  Ny =  5;  % number of elements in 'x'- and 'y'-direction
% hx = 0.1; hy =-0.1; % set explicitly on the code 
% True paraneter values:       Scaled values
% L = 10000 km                 L = 1
% H =  5000 km                 H = 0.5
% E  = 435 GPA                 E  = 1
% nu = 0.2--0.5                nu = 0.2--0.5
% Surf load = 15 MPa           Surf load = ???
% Boxcar = 1000 km             Boxcar = 0.1  (half length)
%
tic
%________________________________ start creating [0,1]^2 domain 
%------------------ simple example: linear in 'y', wh='q0'
%[Node,Edge,Face,Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,...
%      Face_flag,Face_thick,Disco,Discoef,...
%      l_ice,h_ice,rho_ice,rho_earth,...
%      grav,load_surf,...
%      L_char, S_char, U_char, N_char, T_char, Scal,...
%      T_LGM, T_EOG, T,...
%     hx,hy,Nx,Ny,wh] = Square_elast_vq0;
%________________________________ end   creating [0,1]^2 domain 
%________________________________ start creating true glacial domain 
wh = 'gs';
no_domains = 2;
Emagn = 1; % can be 1, 10, 100 (jump in E between the two subdomains)
[L,H,l_ice,h_ice,rho_ice,rho_earth,...
 Disco,Discoef,grav,load_surf,...
 L_char, S_char, U_char, N_char, T_char, Scal,...
 T_LGM, T_EOG, T] = Visco_parameters_new4(no_domains,wh,Emagn); 
% T_LGM, T_EOG, T] = Visco_parameters_Wu(no_domains,wh); 
	  
%[Node,Edge,Face,Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,...
% Face_flag,Face_thick,hx,hy,Nx,Ny] = Rectan_glace(L,H,no_domains); 

%[xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_Wu(L,H);
H0=-max(abs(H));	  
Nx = 6;
Ny = 3;
[xc,yc,hx,hy,Nx,Ny] = Glace_coord_vectors_TH(L,H0,Nx,Ny);	
[Node,Edge,Face,Node_flagx,Node_flagy,...
 Edge_flagx,Edge_flagy,...
 Face_flag,Face_thick] = Rectan_glace_vect(L,H,xc,yc,Nx,Ny,...
                                           no_domains,Disco); 
%________________________________ end   creating true glacial domain 

if(verbose ~= 0)
    figure(1),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,0,16)
end
disp(['Time to create initial mesh: ' num2str(toc)])
% -------------------- Input parameters ---------
levels = input('How many times to refine: ');

   h     = min(abs(hx),abs(hy))/(2^levels);
   sigma = h^2;
%   sigma = 0.1*h^2;  %only for comparisons with Erik
   disp('sigma  = h^2') 

%sigma=h^2/10;
%sigma,wait

solver_type = 1;
   
tic
for lvl=1:levels, 
        [Node,Edge,Face,...
     Node_flagx,Node_flagy,...
     Edge_flagx,Edge_flagy,...
     Face_flag,Face_thick] = my_Refine_quadr(Node,Edge,Face,...
                                   Node_flagx,Node_flagy,...
                                   Edge_flagx,Edge_flagy,...
                                   Face_flag,Face_thick);
%     figure(1),Bvisual_mesh(Node,Edge,Face,1,1,1,0,16)
end

nnodeP = size(Node,2);
nedgeP = size(Edge,2);
nfaceP = size(Face,2);    % number of subdomains (for the pressure)
nallP  = nnodeP;          % number of pressure variables 
disp(['Total number of quadrilaterals for the pressure: ' int2str(nfaceP)])

%figure(1),clf,Bvisual_mesh(Node,Edge,Face,1,1,1,3,16)
disp(['Total number of coarse quadrilaterals: ' int2str(nfaceP)])
%
% ---------------- Create mesh-related matrices -----
%Edge_Node = spalloc(nedgeP,nnodeP,2*nedgeP);
%Face_Edge = spalloc(nfaceP,nedgeP,np*nfaceP);
%for iedge=1:nedgeP,
%    Edge_Node(iedge,Edge(1,iedge))=1;
%    Edge_Node(iedge,Edge(2,iedge))=1;
%end
%for iface=1:nfaceP,
%    Face_Edge(iface,Face(1,iface))=1;
%    Face_Edge(iface,Face(2,iface))=1;	 
%    Face_Edge(iface,Face(3,iface))=1;
%    Face_Edge(iface,Face(4,iface))=1;
%end 
%
%Face_Node_S  = Face_Edge*Edge_Node;  % Face-to-Node-to-Face relation
%for iface=1:nfaceP,
%    rr = find(Face_Node_S(iface,:));
%    Face_Node(1:4,iface,1)=rr';
%nd 
%clear Face_Node_S
%
%Node_Face_no = sum(Face_Node)./2;    % to how many faces each node belongs

% Refine once to obtain the mesh for the displacements

lvl_total = levels + 1;
lvl_coars = max(1,levels-1);
nface_lvl    = zeros(lvl_total,1);
nface_lvl(1) = size(Face,2);
nnode_lvl(1) = size(Node,2);

for lvl=1:1,   % only once
   [Node,Edge,Face,...
    Node_flagx,Node_flagy,...
    Edge_flagx,Edge_flagy,...
    Face_flag,Face_thick,...
    Face_Node,Face_Parent,...
    Face_eorder11,Face_eorder22,...
    nface_lvl,nnode_lvl] = my_Refine_quadr_hier(Node,Edge,Face,...
                                Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,...
                                Face_flag,Face_thick,nface_lvl,nnode_lvl,lvl);
%Bvisual_mesh(Node,Edge,Face,0,0,0,0,16)
end

nnode = size(Node,2);
nedge = size(Edge,2);
nface = size(Face,2);      % number of subdomains (for the displacements)

% ---------------- Create mesh-related matrices (again) -----
Edge_Node = spalloc(nedge,nnode,2*nedge);
Face_Edge = spalloc(nface,nedge,np*nface);
for iedge=1:nedge,
    Edge_Node(iedge,Edge(1,iedge))=1;
    Edge_Node(iedge,Edge(2,iedge))=1;
end
for iface=1:nface,
    Face_Edge(iface,Face(1,iface))=1;
    Face_Edge(iface,Face(2,iface))=1;	 
    Face_Edge(iface,Face(3,iface))=1;
    Face_Edge(iface,Face(4,iface))=1;
end 

% ------------ detect boundary edges
bndry_edge = zeros(nedge,1);
bndry_edge = sum(Face_Edge,1);    % The boundary edges are with sum '1'
[noi,Bndry_Edges]=find(bndry_edge==1); % Bndry_Edges is a list of boundary edges ONLY!
% ------------ detect boundary edges under the ice
if wh=='gs'  % ice for y=0, 0<=x<=length_ice
  [noi,noj]=find(Node(2,:)==0);
  [noii,Bndry_Ice] = find(Node(1,noj)<=l_ice); 
end 

Surface_Nodes = find(Node(2,:)==0); % all surface nodes

Node_Face_noP = sum(Face_Node(:,:,1))./2; %to how many faces each P node belongs
Node_Face_noD = sum(Face_Node(:,:,2))./2; %to how many faces each D node belongs

%depth  = -2e-4;
%[Surface_Nodes, Depth_Nodes] = Level_Nodes(Node,depth);
%Surface_NodesP = Surface_Nodes(1:Nx+(Nx-1)*(2^levels-1));

disp(['End refinement. Elapsed: ' num2str(toc)])
if size(Node,2)<625, figure(1),clf,Bvisual_mesh(Node,Edge,Face,0,0,0,3,12), end

disp(' Ready with the mesh.')

disp(' Start the visco-elastic solver.')
Tmax = T;%input(' Run till time: ');

eps_inner = 5e-1;
Visco_elastic_solverTH_new

return

