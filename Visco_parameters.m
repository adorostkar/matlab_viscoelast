% ----------------------------------------------------------------------
% Parameter setting for the visco-elastic problem
% domains   - number of subdomains (no 1 is the one on the surface ) 
% L         - length (horizontal size) of the domain                (m)
% H         - vertical size of the domain                           (m)
% L_char    - characteristic length of the domain = max(L,H)        (m)
% l_ice     - max width of he ice  sheet                            (m)
% h_ice     - max height of the ice sheet                           (m)
% nju       - Poisson ratio (per subdomain)                dimensionless
% E         - Young modulus (per subdomain)                        (Pa)
% rho_ice   - ice density                                       (kg/m³)
% rho_earth - Earth density                                     (kg/m³)
% eta       - viscosity                                          (Pa s)
% grav      - gravity constant                                   (m/s²)
% S_char    - characteristic stress (S=max(E_i), i = 1:domains)    (Pa)
% U_char    - characteristic displacement                           (m)
% scal_fact - scaling factor (L²/(S*U)                           (m/Pa)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% T_BEG     - time >0, when the ise load starts to be imposed
% T_LGM     - time to reach the glacial maximum
% T_EOG     - time for the ice to melt
% T_run     - time to run after the ice has melted      
% load_surf - load surface = 1(boxcar), 2(ellpce), 3(hyperbola)
%
% The array 'Disco' describes the discontinuity regions
% Disco(*,*,k) is the kth region
% Disco(*,1,k) x-interval 
% Disco(*,2,k) y-interval
% --------------------------
% Discoef(1,k) = nju    in region k
% Discoef(2,k) = E      in region k
% Discoef(3,k) = mu/eta in region k
%
%    --------------------
%    |                  |
%    |   Omega2         | Omega2 - nju and E are variable
%    --------------------
%    |                  |
%    |   Omega1         | Omega2 nju=0.5, E=400e+0 Pa
%    --------------------


%
function [L,H,l_ice,h_ice,rho_ice,rho_earth,...
          Disco,Discoef,grav,load_surf,...
	  L_char, S_char, U_char, N_char, T_char, Scal,...
	  T_LGM, T_EOG, T] = visco_parameters(domains,wh,Emagn)        

E_domains=400e+9*ones(domains,1);
E_domains(2,1)=E_domains(2,1)*Emagn;
% - - - - - - problem and geometry parameters
L0         = 1e+7;   % m  	 
H(1)       =-2e+6;   % m
H(2)       =-4e+6;   % m
l_ice0     = 1e+6;   % m
h_ice0     = 2000;   % m
rho_ice0   =  981;   % kg/m^3
rho_earth0 = 3300;   % kg/m^3
nju(1)     = 0.5;    % dimensionless
nju(2)     = 0.2;    % dimensionless (to be varied)
E0         = max(E_domains); % Pa
mju0       = E_domains./(2*(1+nju'));	 
eta0       = 1.45e+21;  % Pa s
grav       = 9.81;   % m/s^2
spY        = 365*24*60*60;
Years      = 200000; % total simulation time (200000 years)
T0         = Years*spY; % s - total simulation time in seconds

% - - - - - - characteristic values to obtain dimensionless problem
L_char     = abs(L0);  %L_char=max(abs(L),abs(H));
S_char     = E0;       %S_char=max(E) in all subdomains
U_char     = 1; 
N_char     = eta0;% max(of eta) in all subdomains
T_char     = N_char/S_char; % so that S_char*T_char/N_char = 1
Scal       = L_char^2/(S_char*U_char);
%- - - - -
T_LGM0     = 90000*spY;         % Last Glaciation Maximum
T_EOG0     = T_LGM0+10000*spY;  % End Of Glaciation	 
load_surf  = 1; %boxcar


% - - - - - - scaled values
L  = L0/L_char;
H  = H/L_char;
E  = E_domains/S_char;
l_ice     = l_ice0/L_char;
h_ice     = h_ice0/L_char;
rho_ice   = rho_ice0; %??
rho_earth = rho_earth0; %??
eta       = eta0/N_char;
T_LGM     = T_LGM0/T_char;
T_EOG     = T_EOG0/T_char;
T         = T0/T_char;

Disco(1,1,1:domains) =  0;   %form x
Disco(2,1,1:domains) =  L;   %to   x
%  Horizontal split of the domain into strips (two in this case)
Disco(1,2,1) =  0;     %from y
Disco(2,2,1) = H(1);   %to   y
Disco(1,2,2) = H(1);   %from y
Disco(2,2,2) = H(2);   %to   y

Discoef(1,1:domains) = nju(1:domains);
Discoef(2,1:domains) = E';
mju                  = mju0/S_char;
Discoef(3,1:domains) = mju(1:domains)'/eta;% inverse of the Maxwell time 

return
