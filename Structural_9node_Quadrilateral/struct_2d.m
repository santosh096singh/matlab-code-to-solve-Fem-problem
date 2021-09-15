clear all;
E = 30e6; nu = 0.25; strs_flg = 1; 

[C] = constitutive(E, nu, strs_flg);

connect = [1,3,9,7,2,6,8,4,5;
           7,9,15,13,8,12,14,10,11];


% Gausspoints for 3x3 line integration
% --------------------------------
xi = [-0.774597, 0, 0.774597];
xiw = [5/9, 8/9, 5/9];

eta = [-0.774597, 0, 0.774597];
etaw = [5/9, 8/9, 5/9];

% Element 1
coord = [0.0, 0.0;
         3, 0.0;
         3, 1;
         0, 1;
         1.5, 0;
         3,  0.5;
         1.5, 1;
         0, 0.5;
         1.5, 0.5];
     
 % Sum over gauss points:
 ke1 = zeros(18,18);
 for i = 1:3
     for j = 1:3
         keg = elestiff(xi(i), eta(j), coord, C);
         ke1 = ke1 + keg*xiw(i)*etaw(j);
     end
 end
 
 edgeno = 2;
 tvec = [150;0];
 fe1 = zeros(18,1);
 for i = 1:3
    feg = gamat(xi(i), coord, tvec, edgeno);
    fe1 = fe1 + feg*xiw(i);
 end
 
 % Element 2
coord = [0.0, 1.0;
         3, 1.0;
         3, 2;
         0, 2;
         1.5, 1;
         3,  1.5;
         1.5, 2;
         0, 1.5;
         1.5, 1.5];
     
 % Sum over gauss points:
 ke2 = zeros(18,18);
 for i = 1:3
     for j = 1:3
         keg = elestiff(xi(i), eta(j), coord, C);
         ke2 = ke2 + keg*xiw(i)*etaw(j);
     end
 end
 
 edgeno = 2;
 tvec = [100;0];
 fe2 = zeros(18,1);
 for i = 1:3
    feg = gamat(xi(i), coord, tvec, edgeno);
    fe2 = fe2 + feg*xiw(i);
 end
     
% Assembly
K = zeros(30,30);
F = zeros(30,1);

vec = connect(1,:);
for i = 1:9
    for j = 1:9
        K(2*vec(i)-1:2*vec(i), 2*vec(j)-1:2*vec(j)) = K(2*vec(i)-1:2*vec(i), 2*vec(j)-1:2*vec(j)) + ke1(2*i-1:2*i,2*j-1:2*j);
    end
    F(2*vec(i)-1:2*vec(i)) =  F(2*vec(i)-1:2*vec(i)) + fe1(2*i-1:2*i);
end

vec = connect(2,:);
for i = 1:9
    for j = 1:9
        K(2*vec(i)-1:2*vec(i), 2*vec(j)-1:2*vec(j)) = K(2*vec(i)-1:2*vec(i), 2*vec(j)-1:2*vec(j)) + ke2(2*i-1:2*i,2*j-1:2*j);
    end
    F(2*vec(i)-1:2*vec(i)) =  F(2*vec(i)-1:2*vec(i)) + fe2(2*i-1:2*i);
end
        

% Imposition of B.C.
Kreduce = K([3:5,7:24,27:30],[3:5,7:24,27:30]);
Freduce = F([3:5,7:24,27:30]);

% Finding Solution
ureduce = inv(Kreduce)*Freduce;
un = [0;0;ureduce(1:3);0;ureduce(4:21);0;0;ureduce(22:25)];

% Stresses at Gauss Points
Stress_Gauss = zeros(2,9,3);     % No. of elements x No. of gauss points x 3 stress components
for ele = 1:2
    vec = connect(ele,:);
    for i = 1:9
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    ig = 0;
    for i = 1:3
        for j = 1:3
            ig = ig + 1;
            [strs] = stress(xi(i), eta(j), coord, C, uvec');
            Stress_Gauss(ele,ig,1:3) = strs;
        end
    end
end

% Stresses at Nodes
Stress_Nodes = zeros(15,4);    % No. of Nodes x (3 Stress components, no. of occurences)
xi_node = [-1,1,1,-1,0,1,0,-1,0];
eta_node = [-1,-1,1,1,-1,0,1,0,0];
for ele = 1:2
    vec = connect(ele,:);
    for i = 1:9
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    loc_nd = 0;
    for i = 1:9
        nd = vec(i);
        old_stress = Stress_Nodes(nd,1:3);
        ocurrences = Stress_Nodes(nd,4);
        [strs] = stress(xi_node(i), eta_node(i), coord, C, uvec');
        if (ocurrences ~= 0)
            Stress_Nodes(nd,1:3) = (old_stress*ocurrences + strs')/(1+ocurrences);
            Stress_Nodes(nd,4) = ocurrences+1;
        else
            Stress_Nodes(nd,1:3) = strs';
            Stress_Nodes(nd,4) = 1;
        end
    end
end
        

% % % ==============================================
% % %% Printing Intermediate Result to The Output File
% % % ------------------------------------------------
fid=fopen('Steps','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'E = %12.4e, nu = %12.4e, Plane Stress\n\n',E,nu);
fprintf(fid,'Ke1 \n');
fprintf(fid,'----\n\n');
for i = 1:18
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',ke1(i,1:9));
end
fprintf(fid,'\n');
for i = 1:18
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',ke1(i,10:18));
end


fprintf(fid,'Ke2 \n');
fprintf(fid,'----\n\n');
for i = 1:18
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',ke2(i,1:9));
end
fprintf(fid,'\n');
for i = 1:18
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',ke2(i,10:18));
end


fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe1 \n');
fprintf(fid,'-----\n\n');
for i = 1:18
   fprintf(fid,'%10.3e  ',fe1(i));
end
fprintf(fid,'\nfe2 \n');
fprintf(fid,'-----\n\n');
for i = 1:18
   fprintf(fid,'%10.3e  ',fe2(i));
end


fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:30
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',K(i,1:10));
end
fprintf(fid,'\n');
for i = 1:30
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',K(i,11:20));
end
fprintf(fid,'\n');
for i = 1:30
   fprintf(fid,'%10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n\n',K(i,21:30));
end


fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:30
   fprintf(fid,'%12.4e\t',F(i));
end

% fprintf(fid,'\n\nImposition of Boundary Condition\n');
% fprintf(fid,'==================================\n\n');
% fprintf(fid,'K u = F\n');
% fprintf(fid,'--------\n');
% for i = 1:12
%    fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:12),F(i));
% end
% 
fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:25
   fprintf(fid,'%10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e \n\n',Kreduce(i,1:9));
end
fprintf(fid,'\n');
for i = 1:25
   fprintf(fid,'%10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e \n\n',Kreduce(i,10:18));
end
fprintf(fid,'\n');
for i = 1:25
   fprintf(fid,'%10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e   %10.3e\t\t%12.4e\n\n',Kreduce(i,19:25),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'Tn\n');
fprintf(fid,'--\n');
for i = 1:30
   fprintf(fid,'%12.4e\n',un(i));
end

fprintf(fid,'\n\nStresses at Gauss Points\n');
fprintf(fid,'=========================\n\n');
for ele = 1:2
    fprintf(fid,'Element No.: %12d\n',ele);
    for i = 1:4
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    ig = 0;
    for i = 1:3
        for j = 1:3
            ig = ig + 1;
            fprintf(fid,'xi - %10.3e, eta - %10.3e, stresses - %10.3e   %10.3e   %10.3e\n',xi(i),eta(j), Stress_Gauss(ele,ig,1:3));
        end
    end
end
 
fprintf(fid,'\n\nStresses at Nodes\n');
fprintf(fid,'======================\n\n');
for nd = 1:15
    fprintf(fid,'Node:  %12d   Stresses - %10.3e   %10.3e   %10.3e\n',nd, Stress_Nodes(nd,1:3));
end


