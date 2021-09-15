function [C] = constitutive(E, nu, planeflg)

% This function calculate Constitutive Matrix

% INPUT:
% ======
% E = Young's Modulus
% nu = Poisson's ratio
% planeflg = 1  for plane stress
%            2  for plane strain
%
% OUTPUT:
% =======
% C = Constitutive Matrix

if (planeflg == 1) 
   
    C = [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    C = C*E/(1-nu^2);
    
elseif (planeflg == 2)
    
    C = [1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
    C = C*E/(1+nu)/(1-2*nu);
    
end
