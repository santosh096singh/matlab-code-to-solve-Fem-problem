function [fe] = gamat(xi, coord, tvec, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% tvec = traction components
% xi = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge traction is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);

l1 = -0.5*xi*(1-xi);     l2 = 1 - xi^2;     l3 = 0.5*xi*(1+xi);

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = l1; N5 = l2; N2 = l3;
    N3 = 0; N4 = 0; N6 = 0; N7 = 0; N8 = 0; N9 = 0;
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N2 = l1; N6 = l2; N3 = l3;
    N1 = 0; N4 = 0; N5 = 0; N7 = 0; N8 = 0; N9 = 0;    
elseif edgeno == 3
    ledge = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N4 = l1; N7 = l2; N3 = l3;
    N1 = 0; N2 = 0; N5 = 0; N6 = 0; N8 = 0; N9 = 0;
elseif edgeno == 4
    ledge = sqrt((x(4)-x(1))^2 + (y(4)-y(1))^2);
    N1 = l1; N8 = l2; N4 = l3;
    N2 = 0; N3 = 0; N5 = 0; N6 = 0; N7 = 0; N9 = 0;
end

N = zeros(2,18);
N(1,1:2:18) = [N1, N2, N3,N4, N5, N6, N7, N8, N9];
N(2,2:2:18) = [N1, N2, N3,N4, N5, N6, N7, N8, N9];

fe = ledge*N'*tvec*0.5;
  
