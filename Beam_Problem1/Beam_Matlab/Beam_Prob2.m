E = 210e9; I = 0.5*10^-6; q = 2000; P = 10000; M = 3000;

% Element 1
x = [0, 2];
k1 = elestiff(E, I, x);

% Element 2
x = [2, 4];
k2 = elestiff(E, I, x);
f2 = eleload(q,x);

% Element 3
x = [4, 5];
k3 = elestiff(E, I, x);

% Element 4
x = [5, 6];
k4 = elestiff(E, I, x);

% Assembly
K = zeros(10,10);
F = zeros(10,1);

K(1:4,1:4) = k1(1:4,1:4);
K(3:6,3:6) = K(3:6,3:6) + k2(1:4,1:4);
K(5:8,5:8) = K(5:8,5:8) + k3(1:4,1:4);
K(7:10,7:10) = K(7:10,7:10) + k4(1:4,1:4);

F(3:6) = f2(1:4);
F(3) = F(3) + P;
F(8) = F(8) + M;

% Imposition of B.C.
Kreduce = K([2:8,10],[2:8,10]);
Freduce = F([2:8,10]);

% Finding Solution
ureduce = inv(Kreduce)*Freduce;

% Finding Reaction Force
un = [0;ureduce(1:7);0;ureduce(8)];
Fr = K*un;

% FEM displacement
xn = [0,2,4,5,6];
xnume= []; unume = [];
for el = 1:4
    x_n = xn(el:el+1);
    u_n = un((el-1)*2+1:2*(el+1));
    le = x_n(2) - x_n(1);
    xi = [-1:0.2:1]';
    Nx = [(1-xi)/2, (1+xi)/2];
    N1 = (2-3*xi+xi.^3)/4;
    N2 = (1-xi -xi.^2 +xi.^3)/4;
    N3 = (2 + 3*xi -xi.^3)/4;
    N4 = (-1 -xi + xi.^2 + xi.^3)/4;
    Nu = [N1 le*N2/2 N3 le*N4/2];
    xnume = [xnume;Nx*x_n'];
    unume = [unume;Nu*u_n];
end

% Analytical Solution
l = 6; a = 2; b = 4; c = 5; 

xa_1 = [0:0.01:1]*a;
ya_1 = (28.1589*xa_1 - 8.17*xa_1.^3/6)/(E*I);

xa_2 = a + [0:0.01:1]*(b-a);
ya_2 = (28.1589*xa_2 - 8.17*xa_2.^3/6 + 5*(xa_2 - a).^3/3 + (xa_2 - a).^4/12)/(E*I);

xa_3 = b + [0:0.01:1]*(c-b);
ya_3 = (28.1589*xa_3 - 8.17*xa_3.^3/6 + 5*(xa_3 - a).^3/3 + (xa_3 - a).^4/12 - (xa_3 - b).^4/12)/(E*I);

xa_4 = c + [0:0.01:1]*(l-c);
ya_4 = (28.1589*xa_4 - 8.17*xa_4.^3/6 + 5*(xa_4 - a).^3/3 + (xa_4 - a).^4/12 - (xa_4 - b).^4/12 - 3*(xa_4 - c).^2/2)/(E*I);

xana = [xa_1, xa_2, xa_3, xa_4];
yana = [ya_1, ya_2, ya_3, ya_4]*1000;

% Ploting:
h = figure(1);
plot(xana,yana,'b-',xnume,unume,'ro','linewidth',1,'MarkerEdgeColor',...
'k','MarkerFaceColor','r','MarkerSize',8);
legend('Analytical','FEM');
grid on;
set(gcf, 'Position', get(0,'Screensize'));
set(gca,'FontSize',12,'Fontweight','demi');
set(gcf, 'defaultTextInterpreter', 'latex');

% Labelling Axes
xlabel('x (m)','fontsize',18);
ylabel('w (m)','fontsize',18);

% Saving the figure
saveas(h,'beam2','png')

% % ==============================================
% %% Printing Intermediate Result to The Output File
% % ------------------------------------------------

fid=fopen('Steps','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'E = %12.4e, I = %12.4e\n\n',E,I);
fprintf(fid,'Ke1 where xe1 = %12.4e, xe2 = %12.4e, Le = %12.4e\n',xn(1),xn(2),xn(2)-xn(1));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k1(i,1:4));
end
fprintf(fid,'Ke2 where xe1 = %12.4e, xe2 = %12.4e, Le = %12.4e\n',xn(2),xn(3),xn(3)-xn(2));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k2(i,1:4));
end
fprintf(fid,'Ke3 where xe1 = %12.4e, xe2 = %12.4e, Le = %12.4e\n',xn(3),xn(4),xn(4)-xn(3));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k3(i,1:4));
end
fprintf(fid,'Ke4 where xe1 = %12.4e, xe2 = %12.4e, Le = %12.4e\n',xn(4),xn(5),xn(5)-xn(4));
fprintf(fid,'------------------------------------------------------------\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',k4(i,1:4));
end

fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe2 where q = %12.4e, xe1 = %12.4e, xe2 = %12.4e, Le = %12.4e\n',q,xn(2),xn(3),xn(3)-xn(2));
fprintf(fid,'-------------------------------------------------------------------------------\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\n\n',f2(i));
end

fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:10
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:10));
end

fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:10
   fprintf(fid,'%12.4e\n\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:10
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:10),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:8
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',Kreduce(i,1:8),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'wn\n');
fprintf(fid,'--\n');
for i = 1:10
   fprintf(fid,'%12.4e\n\n',un(i));
end

fprintf(fid,'\n\nThe Reaction Forces can be found from\n');
fprintf(fid,'=========================================\n\n');
fprintf(fid,'F = K*wn\n');
fprintf(fid,'---------\n');
for i = 1:10
   fprintf(fid,'%12.4e\n\n',Fr(i));
end
% fprintf(fid,'The Reaction forces are\n');
% fprintf(fid,'=======================\n\n');
% fprintf(fid,'Dof    Reaction force\n');
% fprintf(fid,'===    ==============\n');
% fprintf(fid,'%2i\t\t%10.4e\n',fdof');
% save('filename','fid')
