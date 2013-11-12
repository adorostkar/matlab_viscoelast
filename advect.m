% Advection terms element matrix
% (u \cdot b)\grad v + (\div u) c*v
% A(2*nip,2*nip) in separate displacement ordering
% np is the number of points per f.e.
% k  is the current integration point, to evaluate at

function A = advect(FUN,DER,Coord,np,vec_coeff,nju,wh)

A  =zeros(2*np,2*np);
A11=zeros(np,np); A12=zeros(np,np); A21=zeros(np,np); A22=zeros(np,np);

b1 = b_vec(Coord,1,vec_coeff,wh);
b2 = b_vec(Coord,2,vec_coeff,wh);
if nju==0.5,
   c1 = 0;
   c2 = 0;
else
   c1 = c_vec(Coord,1,vec_coeff,wh);
   c2 = c_vec(Coord,2,vec_coeff,wh);
end

% A11 = -b1*FUN*DER(1,:) + c1*FUN*DER(1,:);
% A12 = -b2*FUN*DER(1,:) + c1*FUN*DER(2,:);
% A21 = -b1*FUN*DER(2,:) + c2*FUN*DER(1,:);
% A22 = -b2*FUN*DER(2,:) + c2*FUN*DER(2,:);

% % To speed up the matlab code, use the fact that b1=c1=0
% % A11 = -b1*FUN*DER(1,:) + c1*FUN*DER(1,:);
% A12 = -b2*FUN*DER(1,:);
% A21 =  c2*FUN*DER(1,:);
% %A22 = -b2*FUN*DER(2,:) + c2*FUN*DER(2,:);
% A22 = FUN*((-b2 + c2)*DER(2,:));

A(1:np,     np+1:2*np) = -b2*FUN*DER(1,:);
A(np+1:2*np,1:np)      =  c2*FUN*DER(1,:);
A(np+1:2*np,np+1:2*np) = FUN*((-b2 + c2)*DER(2,:));
%A = [A11 A12;A21 A22];

return
