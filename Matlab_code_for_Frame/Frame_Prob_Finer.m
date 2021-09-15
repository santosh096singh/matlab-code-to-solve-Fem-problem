E = 70e6; I = 16e-6; A = 8e-3; P = 1000; w = 500;

L1 = 0.2; L2 = 0.5; L3 = 0.08;



% Element 1
x1 = [0, 0, 0, L2/50];
k1 = elestiff(E, A, I, x1);
f1 = eleload(0,w, x1);

% Element 2
x2 = [0, L2/50, 0, L2];
k2 = elestiff(E, A, I, x2);
f2 = eleload(0,w, x2);

% Element 3
x3 = [0, L2, L1/2, L2+L3/2];
k3 = elestiff(E, A, I, x3);

% Element 4
x4 = [L1/2, L2+L3/2, L1, L2+L3];
k4 = elestiff(E, A, I, x4);


% Assembly
K = zeros(15,15);
F = zeros(15,1);

K(1:6,1:6) = k1(1:6,1:6);
K(4:9,4:9) = K(4:9,4:9) + k2(1:6,1:6);
K(7:12,7:12) = K(7:12,7:12) + k3(1:6,1:6);
K(10:15,10:15) = K(10:15,10:15) + k4(1:6,1:6);

F(1:6) = f1(1:6);
F(4:9) = F(4:9)+f2(1:6);
F(8) = F(8) - 2*P;
F(11) = F(11) - 4*P;

% Imposition of B.C.
Kreduce = K(4:12,4:12);
Freduce = F(4:12);

% Finding Solution
ureduce = inv(Kreduce)*Freduce;

% Finding Reaction Force
un = [0;0;0;ureduce;0;0;0];
Fr = K*un;

