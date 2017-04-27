param_chin = [0.0023227, 0.0006304;...
    0.0015850, 	0.0005615];
Kw = 2000; % [kN/m3] Coefficiente di Winkler del terreno
Qk = 35; % [kN/m3] carico dei pallet
L.pl = 5.6; %[m] interasse tra le travi di fondazione (luce libera di inflessione della platea)
L.tr = 5.6; %[m] interasse tra i pali (luce libera di inflessione della trave)
fun = 'palo.Qed(u*1e3) + Kw*u*L.pl*L.tr';   % funzione obiettivo
ris = table;

[rows, ~] = size(param_chin);
for i = 1:rows
    m = param_chin(i,1);
    n = param_chin(i,2);
    palo = chin(m, n);  % legame carico-cedimento alla chin per il palo
    brent('u', 'Fed', [fun ';'], Qk*L.pl*L.tr, 0, 1);
    Rz_palo = palo.Qed(u*1e3);
    K_palo = Rz_palo/u;
    Rz_gr = Kw*u;
    tab_ = table(m, n, u, Rz_palo, Rz_gr, Fed, K_palo, Kw);
    ris = [ris; tab_];
end        
Rz_gr_medio = mean(ris.Rz_gr)
disp(ris)

