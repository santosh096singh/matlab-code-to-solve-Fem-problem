function [ke] = elestiff(xi, eta, coord, C)


% This function calculate element Stiffness Matrix

% INPUT:
% ======
% C = Constitutive Matrix
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix

x = coord(:,1);
y = coord(:,2);

l1 = -0.5*xi*(1-xi);     l2 = 1 - xi^2;     l3 = 0.5*xi*(1+xi);
m1 = -0.5*eta*(1-eta);     m2 = 1 - eta^2;     m3 = 0.5*eta*(1+eta);

l1p = -0.5+xi;    l2p = -2*xi;    l3p = 0.5+xi;
m1p = -0.5+eta;    m2p = -2*eta;    m3p = 0.5+eta;

ddmat = [l1p*m1,   l3p*m1,   l3p*m3,    l1p*m3,  l2p*m1,   l3p*m2,   l2p*m3,  l1p*m2,   l2p*m2;
         l1*m1p,   l3*m1p,   l3*m3p,    l1*m3p,  l2*m1p,   l3*m2p,   l2*m3p,  l1*m2p,   l2*m2p];

J = ddmat*coord;

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = zeros(4,8);
R3(1,1:2:18) = ddmat(1,1:9);
R3(2,1:2:18) = ddmat(2,1:9);
R3(3,2:2:18) = ddmat(1,1:9);
R3(4,2:2:18) = ddmat(2,1:9);

B = R1*R2*R3;
ke = det(J)*(B'*C*B);


