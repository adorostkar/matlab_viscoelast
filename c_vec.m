% -(u \cdot b)\grad v + \div u*c*v
% Definition of the vector field 'c'
% Coord(nip,2)

function c = c_vec(Coord,dim,coeff,wh)

%x = Coord(1);
%Sy = Coord(2);

% glacial rebound
C = coeff(3:4);
c = C(dim);

return
