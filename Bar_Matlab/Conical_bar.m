% Problem Parameters
% ------------------
d1 = 0.02;
d2 = 0.05;
L = 0.8;
rho = 2700;
b = 9.81;
E = 70e9;
g = 9.81;


tht = (d2-d1)/L;
% Gauss Points and weights
% -------------------------
xi1 = -0.774597;   w1 = 5/9;
xi2 = 0.774597;    w2 = 5/9;
xi3 = 0;           w3 = 8/9;

% Element1
% --------
xvec = [0;0.1;0.2]; xvec1 = xvec;
K1 = stiff(xi1,xvec,d1,d2,L,E)*w1 + stiff(xi2,xvec,d1,d2,L,E)*w2 + stiff(xi3,xvec,d1,d2,L,E)*w3;
F1 = loadvec(xi1,xvec,d1,d2,L,rho,b)*w1 + loadvec(xi2,xvec,d1,d2,L,rho,b)*w2 + loadvec(xi3,xvec,d1,d2,L,rho,b)*w3

% Element2
% --------
xvec = [0.2;0.35;0.5];  xvec2 = xvec;
K2 = stiff(xi1,xvec,d1,d2,L,E)*w1 + stiff(xi2,xvec,d1,d2,L,E)*w2 + stiff(xi3,xvec,d1,d2,L,E)*w3
F2 = loadvec(xi1,xvec,d1,d2,L,rho,b)*w1 + loadvec(xi2,xvec,d1,d2,L,rho,b)*w2 + loadvec(xi3,xvec,d1,d2,L,rho,b)*w3

% Element3
% --------
xvec = [0.5;0.65;0.8];  xvec3 = xvec;
K3 = stiff(xi1,xvec,d1,d2,L,E)*w1 + stiff(xi2,xvec,d1,d2,L,E)*w2 + stiff(xi3,xvec,d1,d2,L,E)*w3
F3 = loadvec(xi1,xvec,d1,d2,L,rho,b)*w1 + loadvec(xi2,xvec,d1,d2,L,rho,b)*w2 + loadvec(xi3,xvec,d1,d2,L,rho,b)*w3



% Global Matrices
% ----------------
K = zeros(7,7);
K(1:3,1:3) = K1;
K(3:5,3:5) = K(3:5,3:5) + K2;
K(5:7,5:7) = K(5:7,5:7) + K3;

F = zeros(7,1);
F(1:3) = F1;
F(3:5) = F(3:5) + F2;
F(5:7) = F(5:7) + F3;

% External Load
% -------------
F(3) = F(3) + 20;
F(5) = F(5) - 30;

% Imposition of Boundary Conditions
% ----------------------------------
K_reduce = K(2:7,2:7);
F_reduce = F(2:7);

% Solution
% --------
uvec = inv(K_reduce)*F_reduce;



% Postprocessing
% ---------------
u_node = [0;uvec];
x_node = [0;0.1;0.2;0.35;0.5;0.65;0.8];

% Reaction
% --------
Freac = K*u_node;

xi = [-1:0.25:1]';
N = [-xi.*(1-xi)/2, 1-xi.^2, (xi+1).*xi/2];

xn1 = N*x_node(1:3);
un1 = N*u_node(1:3);
xn2 = N*x_node(3:5);
un2 = N*u_node(3:5);
xn3 = N*x_node(5:7);
un3 = N*u_node(5:7);

xn = [xn1;xn2;xn3];
un = [un1;un2;un3];

%% Analytical Solution
% ====================

delx = 0.0025;

% Deflection in portion 1 due to point load
% -----------------------------------------
P1 = -10;
x = 0:delx:0.2;
up1 = 4*P1*x/pi/E/d1./(d1 + tht*x);
x1 = x;
uA = up1(end);

% Deflection in portion 2 due to point load
% -----------------------------------------
P2 = -30;
x = 0.2:delx:0.5;
xbar = x-0.2;
d1bar = d1 + tht*0.2;
up2 = 4*P2*xbar/pi/E/d1bar./(d1bar + tht*xbar);

x2 = x;
up2 = uA + up2;

uB = up2(end);

% Deflection in element 3 due to point load
% -----------------------------------------
x = 0.5:delx:0.8;
x3 = x;
up3 = uB*ones(1,size(x,2));

% Combining point load deflection for three portions
% --------------------------------------------------
xana = [x1,x2,x3];
up = [up1,up2,up3];

% Analytical deflection due to body force
% ---------------------------------------
dx = d1 + tht*xana;
ug = (rho*g*d2^3/3/E/tht/d1)*xana./dx - rho*g*d1*xana/3/E/tht - rho*g*xana.^2/6/E;

% Total Analytical deflection
% ---------------------------
uana = up + ug;


h = figure(1);
plot(xana,uana,'b-',xn,un,'ro','linewidth',2,'MarkerEdgeColor',...
'k','MarkerFaceColor','r','MarkerSize',8);
hold on;
set(gcf, 'Position', get(0,'Screensize'));
set(gca,'FontSize',12,'Fontweight','demi');
set(gcf, 'defaultTextInterpreter', 'latex');
xlabel('x','fontsize',18);
ylabel('u','fontsize',18);
legend('Analytical','FEM');
grid on
hold on


% Saving the figure
saveas(h,'bar','png')

% % ==============================================
% %% Printing Intermediate Result to The Output File
% % ------------------------------------------------

fid=fopen('Steps','w');
fprintf(fid,'The Element Stiffness matrices are\n','fontweight','bold');
fprintf(fid,'===================================\n');
fprintf(fid,'E = %12.4e\t rho = %12.4e\t b = %12.4e\n\n',E,rho,b);
fprintf(fid,'Ke1 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec1(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',K1(i,1:3));
end
fprintf(fid,'Ke2 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec2(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',K2(i,1:3));
end
fprintf(fid,'Ke3 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec3(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',K3(i,1:3));
end

fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'Fe1 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec1(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',F1(i));
end
fprintf(fid,'Fe2 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec2(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',F2(i));
end
fprintf(fid,'Fe3 where xe = %12.4e,\t%12.4e,\t%12.4e\n',xvec3(1:3));
fprintf(fid,'-------------------------------------------------------\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',F3(i));
end

fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:7));
end

fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\n\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:7
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:7),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:6
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t%12.4e\n\n',K_reduce(i,1:6),F_reduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'un\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\n\n',u_node(i));
end

fprintf(fid,'\n\nThe Reaction Forces can be found from\n');
fprintf(fid,'=========================================\n\n');
fprintf(fid,'F = K*un\n');
fprintf(fid,'---------\n');
for i = 1:7
   fprintf(fid,'%12.4e\n\n',Freac(i));
end