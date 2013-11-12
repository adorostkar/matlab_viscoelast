% --------------------------------------------------------------------
% Preparing matrix to be numerically integrated
% Serves the computation of the element mass matrix (triangles)
%
% FUN   - reference basis functions evaluated at a(n integration) point
% ndof  - total number of degrees of freedom per element (6 in this case)
% nodof - number of degrees of freedom per node (2 in this case)
% --------------------------------------------------------------------


function [MTM]=pre_mass(FUN,ndof,nodof)

MTM = zeros(ndof,ndof);
nod = ndof/nodof;
% % if NOT SDO
% % NT  = zeros(ndof,nodof);
% %       for k=1:nod
% %          for l=1:nodof
% %             NT((k-1)*nodof+l,l) = FUN(k);
% %          end
% %       end
% % 
% % MTM = NT*NT';

% for SDO
M = FUN'*FUN;
MTM = [M zeros(nod); zeros(nod) M];

return