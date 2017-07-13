% dati di resistenza del terreno/pali
clear all
Zona = 'Zona Compattabili';
Caso = 'SLU';
param_chin = [0.0023227, 0.0006304;...
    0.0015850, 	0.0005615];
Kw = 2000; % [kN/m3] Coefficiente di Winkler del terreno

% carichi
carico.pallet = 50; % [kN/m2] carico dei pallet
carico.ppp = 7.5;   % [kN/m2] carico del pacchetto di pavimentazione

% interasse pali
L.pl = 4.35; %[m] interasse tra le travi di fondazione (luce libera di inflessione della platea)
L.tr = 4.84; %[m] interasse tra i pali (luce libera di inflessione della trave)

% dati sezioni in cls
pl = sezione(L.tr*1000, 300);    % sezione della platea
tr = sezione(500, 800 - pl.h);   % sezione della trave ridotta

fun = 'palo.Qed(u*1e3) + Kw*u*L.pl*L.tr';   % funzione obiettivo
ris = table;

% calcolo azioni di progetto
Fd.SLE = (carico.pallet + carico.ppp)*L.pl*L.tr + pl.peso_sezione*L.pl + tr.peso_sezione*L.tr;
Fd.SLU = 1.5*(carico.pallet + carico.ppp)*L.pl*L.tr + 1.3*(pl.peso_sezione*L.pl + tr.peso_sezione*L.tr);

[rows, ~] = size(param_chin);
for i = 1:rows
    m = param_chin(i,1);
    n = param_chin(i,2);
    palo = chin(m, n);  % legame carico-cedimento alla chin per il palo
    brent('u', 'Fed', [fun ';'], Fd.SLU, 0, 1);
    Rz_palo = palo.Qed(u*1e3);
    K_palo = Rz_palo/u;
    Rz_gr = Kw*u;
    Lpl = L.pl;
    Ltr = L.tr;
    tab_ = table(Lpl, Ltr, Fed, m, n, u, Rz_palo, Rz_gr, K_palo, Kw);
    ris = [ris; tab_];
end        
Rz_gr_medio = mean(ris.Rz_gr)
disp(ris)
writetable(ris, [Zona ' - ' Caso '.csv'])
