function ke = elestiff(E, Ie, x)

le = x(2) - x(1);

factor1=Ie*E/le^3;
    ke=[12,6*le,-12,6*le;
        6*le,4*le^2,-6*le,2*le^2;
        -12,-6*le,12,-6*le;
        6*le,2*le^2,-6*le,4*le^2]*factor1;