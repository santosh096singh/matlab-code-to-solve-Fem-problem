function  Ke = stiff(xi,xvec,d1,d2,L,E)

% This function calculate integrand of element stiffness matrix for an 
% element (xvec) at one gauss point (xi)
% =======================================================================
N = [-xi*(1-xi)/2, 1-xi^2, (xi+1)*xi/2];
xe = N*xvec;
tht = (d2-d1)/L;
dx = d1 + tht*xe;
A = pi*dx^2/4;

B = [(xi-0.5), -2*xi, (xi+0.5)];
J = B*xvec;

Ke = A*E*(B'*B)/J;
