k = 1.5;  h = 50; Tinf = 25; Q = 0;

% Gauss Points and weights for four point gauss quadrature
% for triangular area integration
% -------------------------------------------------------
xi1 = 1/3;  eta1 = 1/3;  w1 = -27/48;
xi2 = 0.6;  eta2 = 0.2;  w2 = 25/48;
xi3 = 0.2;  eta3 = 0.6;  w3 = 25/48;
xi4 = 0.2;  eta4 = 0.2;  w4 = 25/48;

% Gausspoints for line integration
% --------------------------------
beta1 = -0.774597;   ew1 = 5/9;
beta2 = 0.774597;    ew2 = 5/9;
beta3 = 0;           ew3 = 8/9;

% Element 1
coord = [0.0, 0.0;
         0.4, 0.0;
         0.4, 0.15;
         0.2, 0.0;
         0.4, 0.075;
         0.2, 0.075];
[ke1dg1,fe1dg1] = domain(xi1, eta1, coord, k, Q);
[ke1dg2,fe1dg2] = domain(xi2, eta2, coord, k, Q);
[ke1dg3,fe1dg3] = domain(xi3, eta3, coord, k, Q);
[ke1dg4,fe1dg4] = domain(xi4, eta4, coord, k, Q);
ke1d = ke1dg1*w1 + ke1dg2*w2 + ke1dg3*w3 + ke1dg4*w4;
fe1d = fe1dg1*w1 + fe1dg2*w2 + fe1dg3*w3 + fe1dg4*w4;

[ke1hg1,fe1hg1] = gamah(beta1, coord, h, Tinf, 2);
[ke1hg2,fe1hg2] = gamah(beta2, coord, h, Tinf, 2);
[ke1hg3,fe1hg3] = gamah(beta3, coord, h, Tinf, 2);
ke1h = ke1hg1*ew1 + ke1hg2*ew2 + ke1hg3*ew3;
fe1h = fe1hg1*ew1 + fe1hg2*ew2 + fe1hg3*ew3;

ke1 = ke1d + ke1h;
%ke1 = ke1d;
fe1 = fe1d + fe1h;


% Element 2
coord = [0.0, 0.0;
         0.4, 0.15;
         0.0, 0.3;
         0.2, 0.075;
         0.2, 0.225;
         0.0, 0.15];
     
[ke2dg1,fe2dg1] = domain(xi1, eta1, coord, k, Q);
[ke2dg2,fe2dg2] = domain(xi2, eta2, coord, k, Q);
[ke2dg3,fe2dg3] = domain(xi3, eta3, coord, k, Q);
[ke2dg4,fe2dg4] = domain(xi4, eta4, coord, k, Q);
ke2d = ke2dg1*w1 + ke2dg2*w2 + ke2dg3*w3 + ke2dg4*w4;
fe2d = fe2dg1*w1 + fe2dg2*w2 + fe2dg3*w3 + fe2dg4*w4;

ke2 = ke2d; fe2 = fe2d;

% Element 3
coord = [0.4, 0.15;
         0.4, 0.3;
         0.0, 0.3;
         0.4, 0.225;
         0.2, 0.3;
         0.2, 0.225];
[ke3dg1,fe3dg1] = domain(xi1, eta1, coord, k, Q);
[ke3dg2,fe3dg2] = domain(xi2, eta2, coord, k, Q);
[ke3dg3,fe3dg3] = domain(xi3, eta3, coord, k, Q);
[ke3dg4,fe3dg4] = domain(xi4, eta4, coord, k, Q);
ke3d = ke3dg1*w1 + ke3dg2*w2 + ke3dg3*w3 + ke3dg4*w4;
fe3d = fe3dg1*w1 + fe3dg2*w2 + fe3dg3*w3 + fe3dg4*w4;

[ke3hg1,fe3hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke3hg2,fe3hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke3hg3,fe3hg3] = gamah(beta3, coord, h, Tinf, 1);
ke3h = ke3hg1*ew1 + ke3hg2*ew2 + ke3hg3*ew3;
fe3h = fe3hg1*ew1 + fe3hg2*ew2 + fe3hg3*ew3;

ke3 = ke3d + ke3h;
%ke3 = ke3d;
fe3 = fe3d + fe3h;


% Assembly
K = zeros(12,12);
F = zeros(12,1);

K([1:3,6,7,11],[1:3,6,7,11]) = ke1(1:6,1:6);
K([1,3,5,11,12,10],[1,3,5,11,12,10]) = K([1,3,5,11,12,10],[1,3,5,11,12,10]) + ke2(1:6,1:6);
K([3,4,5,8,9,12],[3,4,5,8,9,12]) = K([3,4,5,8,9,12],[3,4,5,8,9,12]) + ke3(1:6,1:6);

F([1:3,6,7,11]) = fe1(1:6);
F([1,3,5,11,12,10]) = F([1,3,5,11,12,10]) + fe2(1:6);
F([3,4,5,8,9,12]) = F([3,4,5,8,9,12]) + fe3(1:6);


% Imposition of B.C.
Kreduce = K([1:3,6:8,10:12],[1:3,6:8,10:12]);
Freduce = F([1:3,6:8,10:12]) - (K(4,[1:3,6:8,10:12])*180 + K(5,[1:3,6:8,10:12])*180 + K(9,[1:3,6:8,10:12])*180)';

% Finding Solution
ureduce = inv(Kreduce)*Freduce;
un = [ureduce(1:3);180;180;ureduce(4:6);180;ureduce(7:9)];



