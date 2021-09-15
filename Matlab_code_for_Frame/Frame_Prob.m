E = 70e6; I = 16e-6; A = 8e-3; P = 1000; w = 500;

L1 = 0.2; L2 = 0.5; L3 = 0.08;



% Element 1
x1 = [0, 0, 0, L2];
k1 = elestiff(E, A, I, x1);
f1 = eleload(0,w, x1);

% Element 2
x2 = [0, L2, L1/2, L2+L3/2];
k2 = elestiff(E, A, I, x2);

% Element 3
x3 = [L1/2, L2+L3/2, L1, L2+L3];
k3 = elestiff(E, A, I, x3);


% Assembly
K = zeros(12,12);
F = zeros(12,1);

K(1:6,1:6) = k1(1:6,1:6);
K(4:9,4:9) = K(4:9,4:9) + k2(1:6,1:6);
K(7:12,7:12) = K(7:12,7:12) + k3(1:6,1:6);

F(1:6) = f1(1:6);
F(5) = F(5) - 2*P;
F(8) = F(8) - 4*P;

% Imposition of B.C.
Kreduce = K(4:9,4:9);
Freduce = F(4:9);

% Finding Solution
ureduce = inv(Kreduce)*Freduce;

% Finding Reaction Force
un = [0;0;0;ureduce;0;0;0];
Fr = K*un;

% ==============================================
%% Printing Intermediate Result to The Output File
% ------------------------------------------------

fid=fopen('Steps','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'E = %12.4e, I = %12.4e, A = %12.4e\n\n',E,I,A);
fprintf(fid,'Ke1 where xe1 = %12.4e, ye1 = %12.4e, xe2 = %12.4e, ye2 = %12.4e\n',x1(1:4));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k1(i,1:6));
end
fprintf(fid,'Ke2 where xe1 = %12.4e, ye1 = %12.4e, xe2 = %12.4e, ye2 = %12.4e\n',x2(1:4));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k2(i,1:6));
end
fprintf(fid,'Ke3 where xe1 = %12.4e, ye1 = %12.4e, xe2 = %12.4e, ye2 = %12.4e\n',x3(1:4));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:6
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k3(i,1:6));
end

fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');


for i = 1:6
   fprintf(fid,'%14.4e\n\n',f1(i));
end

fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:12
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:12));
end

fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:12
   fprintf(fid,'%12.4e\n\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:12
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:12),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:6
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',Kreduce(i,1:6),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'wn\n');
fprintf(fid,'--\n');
for i = 1:12
   fprintf(fid,'%12.4e\n\n',un(i));
end

fprintf(fid,'\n\nThe Reaction Forces can be found from\n');
fprintf(fid,'=========================================\n\n');
fprintf(fid,'F = K*wn\n');
fprintf(fid,'---------\n');
for i = 1:12
   fprintf(fid,'%12.4e\n\n',Fr(i));
end
% fprintf(fid,'The Reaction forces are\n');
% fprintf(fid,'=======================\n\n');
% fprintf(fid,'Dof    Reaction force\n');
% fprintf(fid,'===    ==============\n');
% fprintf(fid,'%2i\t\t%10.4e\n',fdof');
% save('filename','fid')
