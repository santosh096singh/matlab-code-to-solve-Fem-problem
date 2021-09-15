function [fe] = gamaq(beta, coord, qn, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% qn = Heat flux 
% beta = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge heat flux is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N = [(1-xi)*(1-2*xi), xi*(2*xi-1), 0, 4*xi*(1-xi), 0, 0];
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N = [0, (1-xi)*(1-2*xi), xi*(2*xi-1), 0, 4*xi*(1-xi),0];
elseif edgeno == 3
    ledge = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
    N = [xi*(2*xi-1), 0, (1-xi)*(1-2*xi), 0, 0, 4*xi*(1-xi)];
end

fe = N'*qn*ledge*0.5;
  
