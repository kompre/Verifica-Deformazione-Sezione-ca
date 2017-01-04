clearvars
clc
%% Dimensionamento Platea
% in condizioni SLE la platea � soggeta al carico q = 66.3 kN/m^2. Il
% criterio dimensionante � la freccia nel punto pi� sfavorevole, i.e. nel
% punto medio della platea compreso tra due pali. Pertanto la luce
% effettiva di calcolo � pari a (Lpl + Ltr)^2.

% freccia limite
f_lim = 1/3000;

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

%% LIMITAZIONI AL PROBLEMA DI OTTIMIZZAZIONE
% Di seguito si elencano le limitazioni e le costanti che caratterizzano
% questo problema di ottimizzazione. Sono inserite le dimensioni della
% fondazione, i limiti per il numero di travi (i.e. luce libera di
% inflessione della platea, i limiti per il numero di pali/trave (i.e. luce
% libera di inflessione delle travi).
% In questa parte si indicano i costi unitari al peso ed al palo.
L.X = 56;   % [m] dimensione in x della fondazione
L.Y = 79;   % [m] dimensione in y della fondazione
num.tr = 10:30; % numeri di travi nella fondazione
num.pali_tr = 14:40; % numero di pali per trave
cu.steel = 0.80; % costo unitatio acciaio [�/kg]
cu.palo  = 600;  % costo unitario palo [�/kg]
peso_acciaio = 7850;    % [kg/m^3] peso unitario dell'acciaio

%%
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

%% Funzione freccia elastica in mezzeria
f_el = @(q,l,E,J) -1/384 * q*l^4/(E*J);

%% Variazione luce libera di inflessione della platea
L.pl = L.X ./ (num.tr - 1);      % [m] variazione della luce libera di inflessione della platea
L.tr = L.Y ./ (num.pali_tr - 1); % [m] variazione della luce libera di inflessione della trave

%% 
num.pali = num.tr' * num.pali_tr;   % numero totale di pali nelle varie configurazioni size: lenght(num.tr) x length(numpali_tr)



%% Inizio ciclo principale

% preallocamento variabili
arm.pl = repmat(arm.pl, length(L.pl),1);   % prealloco la struttura arm.pl (N)
arm.tr = repmat(arm.tr, length(L.pl), length(L.tr)); % prealloco la struttura arm.tr  (NxM)
M_max.pl = zeros(size(L.pl));
M_max.tr = zeros(length(L.pl), length(L.tr));
ris.pl.s_min = zeros(size(length(L.pl)));
ris.tr.s_min = zeros(length(L.pl), length(L.tr));
soll.tr = repmat(soll.tr, length(L.pl), 1);

for i_pl = 1:length(L.pl)
    % dimensionamento platea
    M_max.pl(i_pl) = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
    As_min = M_max.pl(i_pl)*1e6 / (0.9*max(arm.pl(i_pl).d) * mat.steel.fyd) /arm.pl(i_pl).nb(2);   % armatura minima con criterio semplificato dimensionata allo SLU
    diam_min = ceil( ceil( sqrt( 4*As_min/pi ) ) /2 ) * 2;   % arrotonda ad un diametro reale
    arm.pl(i_pl).diam = [12;diam_min];
    % calcolo freccia
    ris.pl.s_min(i_pl) = calcoloFreccia(L.pl(i_pl), geom.pl, arm.pl(i_pl), mat, soll.pl, cost, dx, 'semplificato', f_el);
    
    % aggiornamento variabile di carico agente sulla trave
    soll.tr(i_pl).q = soll.pl.q * L.pl(i_pl);
    soll.tr(i_pl).qSLU = soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
    
    % dimensionamento della trave
    for i_tr = 1:length(L.tr)        
        M_max.tr(i_pl, i_tr) = soll.tr(i_pl).M( soll.tr(i_pl).qSLU, L.tr(i_tr), L.tr(i_tr)/2);
        As_min = M_max.tr(i_pl, i_tr)*1e6 / (0.9*max(arm.tr(i_pl, i_tr).d) * mat.steel.fyd) /arm.tr(i_pl, i_tr).nb(2);   % armatura minima con criterio semplificato dimensionata allo SLU
        diam_min = ceil( ceil( sqrt( 4*As_min/pi ) ) /2 ) * 2;   % arrotonda ad un diametro reale
        arm.tr(i_pl, i_tr).diam = [diam_min;14;14;diam_min];
        % calcolo freccia della trave
        ris.tr.s_min(i_pl, i_tr) = calcoloFreccia(L.tr(i_tr), geom.tr, arm.tr(i_pl, i_tr), mat, soll.tr(i_pl), cost, dx, 'semplificato', f_el);
        
    end
end
%% Combinazione dei risultati
for i_pl = 1:length(L.pl)
    for i_tr = 1:length(L.tr)
        L.tot(i_pl, i_tr) = sqrt(L.pl(i_pl)^2 + L.tr(i_tr)^2);  % lunghezza combinata su cui verificare la freccia limite
        ris.tot.s_min(i_pl, i_tr) = ris.pl.s_min(i_pl) + ris.tr.s_min(i_pl, i_tr);
        ris.tot.s_rel(i_pl, i_tr) = ris.tot.s_min(i_pl, i_tr)/(L.tot(i_pl, i_tr)*1e3);    % freccia in termini di luce libera di inflessione      
    end
end
