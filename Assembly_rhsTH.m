% Given: Node,Edge
% asembling of the right-hand-side vector (surface load only)
% The load is time-dependent

function [rhs_d,rhs_p] = Assembly_rhsTH(Node,Edge,...
                                      wh,l_ice,h_ice,rho_ice,grav,Bndry_Ice,...
                                      L_char,S_char,U_char,nnodeP,...
                                      time_length,T_BEG,T_LGM,T_EOG)

nedge = size(Edge,2);   % number of edges
nnode = size(Node,2);   % number of nodes

rhs_d = zeros(2*nnode,1);
rhs_p = zeros(nnodeP,1);
% - - - - - - Add surface forces to the right-hand side
if ~(strcmp(wh,'g1')&strcmp(wh,'g2')), % no surface forces,exact sol known
% 
  for ke= 1:length(Bndry_Ice),
      k = Bndry_Ice(ke);
       [G0,len] = lin_num_int2L_var(Node(1,Edge(1,k)),Node(2,Edge(1,k)),...
                      Node(1,Edge(2,k)),Node(2,Edge(2,k)),...
				     wh,l_ice,h_ice,rho_ice,grav,...
				     time_length,T_BEG,T_LGM,T_EOG);
       G = (L_char/S_char)*G0;   
       rhs_d(Edge(1,k),1)       = rhs_d(Edge(1,k),1)       + G(1);
       rhs_d(Edge(2,k),1)       = rhs_d(Edge(2,k),1)       + G(1);
       rhs_d(Edge(1,k)+nnode,1) = rhs_d(Edge(1,k)+nnode,1) + G(2);
       rhs_d(Edge(2,k)+nnode,1) = rhs_d(Edge(2,k)+nnode,1) + G(2);

    end
end
return
