% -(u \cdot b)\grad v + (\div u) c*v
% Definition of vector field 'b'
% Coord(nip,2)

function b = b_vec(Coord,dim,coeff,wh)

%x = Coord(1);
%y = Coord(2);

% glacial rebound
B = coeff(1:2);
b = B(dim);

return
