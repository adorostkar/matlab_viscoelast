% ------------------ symbolic form
% Quadrilateral elements, bilinear b.f. 
% Modified Taylor-Hood elements
% Coord(4,2)

function [B1_elem,B2_elem]=TH_symbolic(CoordP,Coord)

x  = sym('x');
y  = sym('y');
I4 = eye(4);
    
for k=1:4
    QP(k,1:4) = [CoordP(k,1) CoordP(k,2) CoordP(k,1)*CoordP(k,2) 1];
end
cP = QP\I4; % the columns contain the coefficients of the coarce b.f.

for k=1:4,
    fiP(k) = cP(1,k)*x+cP(2,k)*y+cP(3,k)*x*y+cP(4,k);
    fiPdx(k) = diff(fiP(k),'x');
    fiPdy(k) = diff(fiP(k),'y');
end


for k=1:4
    Q(k,1:4) = [Coord(k,1) Coord(k,2) Coord(k,1)*Coord(k,2) 1];
end
c = Q\I4; % the columns contain the coefficients of the coarce b.f.

for k=1:4,
    fi(k) = c(1,k)*x+c(2,k)*y+c(3,k)*x*y+c(4,k);
    fidx(k) = diff(fi(k),'x');
    fidy(k) = diff(fi(k),'y');
end
fx = min(Coord(:,1)); tx = max(Coord(:,1));
fy = min(Coord(:,2)); ty = max(Coord(:,2));
for k=1:4
    for l=1:4
        B1sym(l,k) = int(int(fidx(l)*fiP(k),y,fy,ty),x,fx,tx);
        B2sym(l,k) = int(int(fidy(l)*fiP(k),y,fy,ty),x,fx,tx);
    end
end   
double(B1sym)
double(B2sym)

% Laplacian part
for k=1:4
    for l=1:4
        Lsym(l,k) = int(int((fidx(l)*fidx(k)+fidy(l)*fidy(k)),y,fy,ty),x,fx,tx);
    end
end   

% Rot part
for k=1:4
    for l=1:4
        R11sym(l,k) = -int(int(fidy(l)*fidy(k),y,fy,ty),x,fx,tx);
	R12sym(l,k) =  int(int(fidx(l)*fidy(k),y,fy,ty),x,fx,tx);
        R21sym(l,k) =  int(int(fidy(l)*fidx(k),y,fy,ty),x,fx,tx);
	R22sym(l,k) = -int(int(fidx(l)*fidx(k),y,fy,ty),x,fx,tx);
    end
end 

Rsym= [R11sym R12sym;R21sym R22sym];


