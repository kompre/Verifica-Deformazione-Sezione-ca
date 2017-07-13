f_el = @(q,l,E,J) -1/384 * q*l^4/(E*J);
dx = 0.001;
arm_.nb =[4;2;2;4;3]
arm_.diam = [22;14;14;18;22]
arm_.d = [250;280;520;750;750]
f = calcoloFreccia(L.pl, geom.pl, arm_, fck, cost.c, dx, soll.pl, 'semplificato', f_el);