function fe = eleload(q, x)

le = x(2) - x(1);
fe =  [q*le/2 ; q*le^2/12 ; q*le/2 ; -q*le^2/12];