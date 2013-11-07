% Gives the value of the body force
% wh is a problem identifier
% Important: F is in 'separate displacement ordering' !!
% ---------------------------------------------------------
% ---------------------------------------------------------
% Dimensions: the domain   - m
%             gravity      - N/m^2
%     specific weight (ro) - kg/m**3 (otnositelno teglo)
%             body_force   - 
%
% Young's module (E): beton: 30.0*10^9 N/m^2 
%                     rock:   6.0*10^9 N/m^2
% Poisson ratio (nu): beton: 0.2            
%                     rock:  0.25
% Surface force: p=ro_water*g*h N/m^2
%                g=9.80665 N/m^2 (gravity constant)
%     water dencity=1.0*10^3 kg/m^3 (specific weight of the water)     
%     beton dencity=2.2*10^3 kg/m^3 (specific weight of the beton) 
%     rock  density=.....    
% ---------------------------------------------------------

function [F,FP1] = body_force_sdo(face_flag,thickness,Coord,...
                                  mju,rho,coeff,wh)

F = zeros(8,1); FP1 = zeros(4,1);
if (wh~='g1')&(wh~='D2'),
  return % no body forces
else 
  if wh=='g1',     %computes \int F*phi and \int nabla_F*phi
     for k=1:4
        F(k,1)  =       g1_force(Coord(k,:),mju,rho,H,L,q,r,1,coeff,wh);
        F(k+3,1)=       g1_force(Coord(k,:),mju,rho,H,L,q,r,2,coeff,wh);
	FP1(k,1)= g1_nabla_force(Coord(k,:),mju,rho,H,L,q,r,  coeff,wh);      
     end
   
  end 
  if wh=='D2',
     gravity = 9.80665; 
     if face_flag==2,
        ro = 2.4*10^(-3)*thickness; %kg/m^3 beton 
     else
        ro = 2.5*10^(-3)*thickness; %kg/m^3 rock foundation
     end
     F(4,1) = -gravity*ro;
     F(5,1) = -gravity*ro;
     F(6,1) = -gravity*ro;
  end
end

return
