function fe = eleload(q_axial,q_transverse, x)

% This function calculate load vector on a frame element for a uniform q

lx = x(3) - x(1);
ly = x(4) - x(2);
le = sqrt(lx^2 + ly^2);

l = lx/le;  m = ly/le; lm = l*m;


fe_beam =  [q_transverse*le/2 ; q_transverse*le^2/12 ; q_transverse*le/2 ; -q_transverse*le^2/12];
fe_bar = (q_axial*le)*[0.5;0.5];

fep = zeros(6,1);
fep([1;4]) = fe_bar;
fep([2;3;5;6]) = fe_beam;

Q = zeros(6,6);
Q(1:2,1:2) = [l,m;m,-l];
Q(4:5,4:5) = [l,m;m,-l];
Q(3,3) = 1;
Q(6,6) = 1;

fe = Q'*fep;