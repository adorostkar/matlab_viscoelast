% Initial quadrilateral mesh
% b.c. such that the solution is displ_x=0; displ_y - linear
%
% Face_flag(iface)=disco_region
%---->  Regular meshsize in each direction
% Data: Nx - number of intervals in x-direction
%       Ny - number of intervals in y-direction
%       [hx] - variable stepsize in x-direction
%       [hy] - variable stepsize in y-direction

function [Node,Edge,Face,...
          Node_flagx,Node_flagy,Edge_flagx,Edge_flagy,...
         Face_flag,Face_thick] = Rectan_glace_vect(L,H,xc,yc,Nx,Ny,...
                                                   domains,Disco,wh)
% ------- subdomain 1:        
%   Node
nnode = Nx*Ny;
Node  = zeros(2,nnode);
for l=1:Ny,
    for k=1:Nx,
        Node(1,(l-1)*Nx + k) = xc(k);
        Node(2,(l-1)*Nx + k) = yc(l);
    end
end

% Edge
nedge = Nx*(Ny-1)+Ny*(Nx-1);
Edge  = zeros(2,nedge);
edge_no=0;
for l=1:Ny-1,
    for k=1:Nx-1
        edge_no = edge_no + 1;
        Edge(1,edge_no) = (l-1)*Nx + k;  
        Edge(2,edge_no) = (l-1)*Nx + k + 1;
    end
    for k=1:Nx  
        edge_no = edge_no + 1;
        Edge(1,edge_no) = (l-1)*Nx + k;    
        Edge(2,edge_no) = l*Nx     + k;
    end
end
l = Ny;
for k=1:Nx-1
    edge_no = edge_no + 1;
    Edge(1,edge_no) =  (l-1)*Nx + k;       
    Edge(2,edge_no) =  (l-1)*Nx + k + 1;
end

% Face
nface = (Nx-1)*(Ny-1);
Face  = zeros(4,nface);
for l=1:Ny-1,
    for k=1:Nx-1
        face_no = (l-1)*(Nx-1) + k;
        iedge   = (l-1)*(2*Nx-1);
        Face(1,face_no) = iedge            + k;
        Face(2,face_no) = iedge + (Nx-1)   + k;
        Face(3,face_no) = iedge + (2*Nx-1) + k;
        Face(4,face_no) = iedge + (Nx-1)   + k + 1;
        all_y_coord = Node(2,Edge(:,Face(:,face_no)));
        if min(abs(all_y_coord))>=abs(Disco(2,2,1)),
           Face_flag(face_no,1)  = 2;
        else
           Face_flag(face_no,1)  = 1;
        end   
    end
end    

%Bvisual_mesh(Node,Edge,Face,1,1,1,3,16)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nnode = size(Node,2);
Node_flagx(1:nnode,1)  = 0;
Node_flagy(1:nnode,1)  = 0;
if wh == 'g0', % benchmark elasticity
   G1=find(Node(1,:)==0);
   G2=find(Node(1,:)==1);
   G3=find(Node(2,:)==0);
   Gy=find(Node(2,:)==-1);
   Gx=cat(2,G1,G2,Gy);
   Node_flagx(Gx,1) = Flags('Dirichlet');
   Node_flagy(Gy,1) = Flags('Dirichlet');   
else
   for k=1:nnode,
       if (Node(1,k)==0)|(Node(2,k)==H(domains)),%|(Node(1,k)==L),
          Node_flagx(k,1)  = Flags('Dirichlet');
       end
       if (Node(2,k)==H(domains)),%|(Node(1,k)==L),
          Node_flagy(k,1)  = Flags('Dirichlet');
       end
   end
end
   nedge = size(Edge,2);
   Edge_flagx(1:nedge,1) = 0;
   Edge_flagy(1:nedge,1) = 0;
   
if wh == 'g0', % benchmark elasticity
   for k=1:nedge,
        k1 = Edge(1,k); k2=Edge(2,k);
       if (ismember(k1,Gx))&(ismember(k2,Gx)),
              Edge_flagx(k,1)  = Flags('Dirichlet');
       end
       if (ismember(k1,Gy))&(ismember(k2,Gy)),
              Edge_flagy(k,1)  = Flags('Dirichlet');
       end 
   end
else
   for k=1:nedge,
       k1 = Edge(1,k); k2=Edge(2,k);
       if ((Node(2,k1)==H(domains))&(Node(2,k2)==H(domains)))|...
          ((Node(1,k1)==0)&(Node(1,k2)==0))|...
          ((Node(1,k1)==L)&(Node(1,k2)==L)),
              Edge_flagx(k,1)  = Flags('Dirichlet');
       end
       if ((Node(2,k1)==H(domains))&(Node(2,k2)==H(domains)))|...
          ((Node(1,k1)==L)&(Node(1,k2)==l)),
              Edge_flagy(k,1)  = Flags('Dirichlet');
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nface = size(Face,2);
Face_thick(1:nface,1) = 1;

return
