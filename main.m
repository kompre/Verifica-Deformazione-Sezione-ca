clearvars
clc
%% Dimensionamento Platea
% in condizioni SLE la platea è soggeta al carico q = 66.3 kN/m^2. Il
% criterio dimensionante è la freccia nel punto più sfavorevole, i.e. nel
% punto medio della platea compreso tra due pali. Pertanto la luce
% effettiva di calcolo è pari a (Lpl + Ltr)^2.

% Dati geometrici sezione della platea
geom.pl.b = 1000;
geom.pl.h = 300;
geom.pl.x0 = 0;
geom.pl.y0 = 0;
% Dati geometrici sezione della trave
geom.tr.b = 500;
geom.tr.h = 800;
geom.tr.x0 = 0;
geom.tr.y0 = 0;

% dati materiale
mat.cls = derivaCaratteristicheCA(25,30);
mat.steel = derivaCaratteristicheAcciaio;

% dati sollecitazione della platea
soll.pl.q = 66.3;  % [kN/m2] carico agente in condizioni SLE
soll.pl.qSLU = 96.19;  % [kN/m2] carico agente in condizioni SLU
soll.pl.N = 0;
soll.pl.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata
% dati sollecitazione della trave
soll.tr.q = 0;  % inizializzo il valore di carico da stimare per ogni iteraizione in funzione della luce di inflessione della platea
soll.tr.qSLU = 0;  % idem come sopra
soll.tr.N = 0;
soll.tr.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata

% costanti
cost.c = 0.5;
cost.rot = 0;
cost.spo = 0;
% passo di integrazione
dx = 0.001;
% dati armatura platea tipo
arm.pl.nb = [5;5];
arm.pl.diam = [12;12];
arm.pl.d = [40;260];
% dati armatura trave tipo
arm.tr.nb = [4;2;2;4];
arm.tr.diam = [14;14;14;14];
arm.tr.d = [50;270;530;750];

%% VARIAZIONE DELLA LUCE LIBERA DI INFLESSIONE DELLA PLATEA
L.pl = 1:0.1:6; % [m] variazione della luce libera di inflessione della platea con risoluzione 10 cm
%% VARIAZIONE DELLA LUCE LIBERA DI INFLESSIONE DELLA TRAVE
L.tr = 1:0.1:10; % [m] variazione della luce libera di inflessione della platea con risoluzione 10 cm

arm.pl = repmat(arm.pl, length(L.pl),1);   % prealloco la struttura arm.pl (N)
arm.tr = repmat(arm.tr, length(L.pl), length(L.tr)); % prealloco la struttura arm.tr  (NxM)

for i_pl = 1:length(L.pl)
    tic
    % dimensionamento platea
    M_max.pl(i_pl) = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
    As_min = M_max.pl(i_pl)*1e6 / (0.9*max(arm.pl(i_pl).d) * mat.steel.fyd / mat.steel.gamma_s) /arm.pl(i_pl).nb(2);   % armatura minima con criterio semplificato dimensionata allo SLU
    diam_min = ceil( ceil( sqrt( 4*As_min/pi ) ) /2 ) * 2;   % arrotonda ad un diametro reale
    arm.pl(i_pl).diam = [12;diam_min];
    % calcolo freccia
    ris.pl.s_min(i_pl) = calcoloFreccia(L.pl(i_pl), geom.pl, arm.pl(i_pl), mat, soll.pl, cost, dx);
    
    % dimensionamento trave

    soll.tr.q = soll.pl.q * L.pl(i_pl);
    soll.tr.qSLU = soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
    for i_tr = 1:length(L.tr)        
        M_max.tr(i_pl,i_tr) = soll.tr.M(soll.tr.qSLU, L.tr(i_tr), L.tr(i_tr)/2);
        As_min = M_max.tr(i_pl, i_tr)*1e6 / (0.9*max(arm.tr(i_pl, i_tr).d) * mat.steel.fyd / mat.steel.gamma_s) /arm.tr(i_pl, i_tr).nb(2);   % armatura minima con criterio semplificato dimensionata allo SLU
        diam_min = ceil( ceil( sqrt( 4*As_min/pi ) ) /2 ) * 2;   % arrotonda ad un diametro reale
        arm.tr(i_pl, i_tr).diam = [diam_min;14;14;diam_min];
        % calcolo freccia
        ris.tr.s_min(i_pl, i_tr) = calcoloFreccia(L.tr(i_tr), geom.tr, arm.tr(i_pl, i_tr), mat, soll.tr, cost, dx);

    end
    toc
end
%% Stima dell'armatura
% l'armatura è stimata con metodo semplificato in funzione del massimo
% momento in campata






%%

