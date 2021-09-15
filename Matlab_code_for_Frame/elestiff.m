function ke = elestiff(E, A, Ie, x)

lx = x(3) - x(1);
ly = x(4) - x(2);
le = sqrt(lx^2 + ly^2);

l = lx/le;  m = ly/le; lm = l*m;

factor1=Ie*E/le^3;


ke_beam=[12,6*le,-12,6*le;
    6*le,4*le^2,-6*le,2*le^2;
    -12,-6*le,12,-6*le;
    6*le,2*le^2,-6*le,4*le^2]*factor1;

ke_bar = (E*A/le)*[1,-1;-1,1];

kep = zeros(6,6);
kep(1:3:4,1:3:4) = ke_bar;
kep([2:3,5:6],[2:3,5:6]) = ke_beam;

Q = zeros(6,6);
Q(1:2,1:2) = [l,m;m,-l];
Q(4:5,4:5) = [l,m;m,-l];
Q(3,3) = 1;
Q(6,6) = 1;

ke = Q'*kep*Q;
    
       
    
    
    