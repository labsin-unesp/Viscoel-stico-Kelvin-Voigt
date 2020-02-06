function [u,ve,ac] = newmark(K,C,M,f,dt)

% Newmark time-integration scheme
% ----------------------------------------------------
%
% PURPOSE: Numerical integration of
%   M*a(t) + C*v(t)  + K*u(t) = f(t)
% 	
% INPUT:   
%   K  : system stiffness matrix   
%	C  : system damping matrix 
%	M  : system mass matrix 
%	f  : system load vector
%	dt : time integration step
%           
% OUTPUT:  
%	ac : accleration at nodal variables 
%	ve : velocity at nodal variables 
%	u  : displacement at nodal variables 
%
% ----------------------------------------------------
% 27.11.2002 Lars Andersen, AAU

% Set the integration parameters ---------------------

gamma = 1/2 ; beta = 1/4 ;

% Number of time steps -------------------------------

n = size(f,2) ;

% Initialize the vectors u, v and a ------------------

u = zeros(size(f)) ;
ve = zeros(size(f)) ;
ac = zeros(size(f)) ;

% Calculate the inverse system matrix ----------------

M1 = inv(M + gamma*C*dt + beta*K*dt^2) ;

% Perform loop over the time steps

for j = 1 : n-1
   
   % Predicted values of u and v ----------------------
   
   vv = ve(:,j) + ac(:,j)*dt ;
   uu = u(:,j) + ve(:,j)*dt + 0.5*ac(:,j)*dt^2 ;
   
   % Corrected values ---------------------------------
   
   ac(:,j+1) = ac(:,j) + M1*(f(:,j+1)-M*ac(:,j)-C*vv-K*uu); 
   ve(:,j+1) = vv + gamma*(ac(:,j+1)-ac(:,j))*dt;
   u(:,j+1) = uu + beta*(ac(:,j+1)-ac(:,j))*dt^2;
   
end

% End of file ----------------------------------------