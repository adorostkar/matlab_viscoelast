% Advection terms element matrix
% (u \cdot b)\grad v + (\div u) c*v
% A(2*nip,2*nip) in separate displacement ordering
% np is the number of points per f.e.
% k  is the current integration point, to evaluate at

function A = advect1(FUN,DER,Coord,np,coeff,nju,wh)

%A = zeros(2*np,2*np);
A11=zeros(np,np); A12=zeros(np,np); A21=zeros(np,np); A22=zeros(np,np);

b1 = b_vec(Coord,1,coeff,wh);
b2 = b_vec(Coord,2,coeff,wh);
if nju==0.5,
   c1 = 0;
   c2 = 0;
else
   c1 = c_vec(Coord,1,coeff,wh);
   c2 = c_vec(Coord,2,coeff,wh);
end

A11 = -b1*DER(1,:)'*FUN + c1*DER(1,:)'*FUN;
A12 = -b2*DER(1,:)'*FUN + c1*DER(2,:)'*FUN;
A21 = -b1*DER(2,:)'*FUN + c2*DER(1,:)'*FUN;
A21 = -b2*DER(2,:)'*FUN + c2*DER(2,:)'*FUN;

A = [A11 A12;A21 A22];

return
