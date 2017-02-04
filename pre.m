clearvars
close all
clc
%% Dimensionamento Platea
% in condizioni SLE la platea è soggeta al carico q = 66.3 kN/m^2. Il
% criterio dimensionante è la freccia nel punto più sfavorevole, i.e. nel
% punto medio della platea compreso tra due pali. Pertanto la luce
% effettiva di calcolo è pari a (Lpl + Ltr)^2.

% freccia limite
f_lim = 1/2500;

% Dati geometrici sezione della platea
geom.pl = sezione(1000, 300, 0, 0, 1, 100);

% Dati geometrici sezione della trave
geom.tr = sezione(500, 800, 0, 0, 1, 100);


% dati materiale
fck = 25;

% dati sollecitazione della platea
% composizione del carico:
%   G2 = 8.8 kN/m2 
%   Qk = 50.0 kN/m2
sollecitazioni.pl.q = 50;  % [kN/m2] carico agente in condizioni SLE
sollecitazioni.pl.qSLU = 1.5*50 + 1.3*8.8;  % [kN/m2] carico agente in condizioni SLU
sollecitazioni.pl.N = 0;
sollecitazioni.pl.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata
sollecitazioni.pl.V = @(q,l,x) -q*x + q*l/2; % taglio lungo l'asse della trave
% dati sollecitazione della trave
sollecitazioni.tr.q = 0;  % inizializzo il valore di carico da stimare per ogni iterazione in funzione della luce di inflessione della platea
sollecitazioni.tr.qSLU = 0;  % idem come sopra
sollecitazioni.tr.N = 0;
sollecitazioni.tr.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata
sollecitazioni.tr.V = @(q,l,x) -q*x + q*l/2;


%% LIMITAZIONI AL PROBLEMA DI OTTIMIZZAZIONE
% Di seguito si elencano le limitazioni e le costanti che caratterizzano
% questo problema di ottimizzazione. Sono inserite le dimensioni della
% fondazione, i limiti per il numero di travi (i.e. luce libera di
% inflessione della platea, i limiti per il numero di pali/trave (i.e. luce
% libera di inflessione delle travi).
% In questa parte si indicano i costi unitari al peso ed al palo.
L.X = 56.46;   % [m] dimensione in x della fondazione
L.Y = 79.44;   % [m] dimensione in y della fondazione
L.b_max = 12; % [m] lunghezza massima delle barre

% La platea in oggetto è caratterizzata da 3 campi definiti dai plinti dei
% pilastri del magazzino che costituiscono una invariante del problema. Le
% dimensioni di questi campi sono 17.35, 21.45, 16.80, pertanto si è deciso
% che i campi di larghezza 17.35 e 16.80 abbiano lo stesso numero di
% campate, mentre il campo di larghezza 21.45 avrà una campata in più
% (relazione valida per un range limitato di variazione delle campate
num.campate = 3*(3:6)+1;
num.pali_tr = 10:20; % numero di pali per trave 14:40

% costi unitati
cu.cls = 83;    % [€/m3] costo del cls
cu.pompaggio = 3.64;    % [€/m3] costo del pompaggio del cls
cu.casseforma = 15; % [€/m2] costo dei casseri (da applicare alle travi)
cu.steel = 0.79; % costo unitatio acciaio [€/kg]
cu.palo  = 670.08;  % costo unitario palo [€/kg]

peso_acciaio = 7850;    % [kg/m^3] peso unitario dell'acciaio

RZ.palo.max = 1150; % [kN] carico massimo ammissibile per il palo
Kw.ground = 2000; % [kN/m3] coefficiente di Winkler del terreno
Kw.palo = 125e3;    % [kN/m] modulo secante per legame palo terreno

%%
% costanti
cost.k = 4*sqrt(6)/9; % amplificatzione della freccia per ottenerre il massimo rapporto df/dL in funzione delle condizioni di appoggio (1 per app-app, 4*sqrt(6)/9 per inc-inc)
f_lim = f_lim/cost.k;
cost.c = 0.5;   % costante per il calcolo del parametro zeta (c = 0.5 per carichi ripetuti o lunga durata, c = 1 per carichi brevi)
cost.rot = 0;   % costante di rotazione (integrale della curvatura)
cost.spo = 0;   % costante di spostamento (integrale della rotazione)

% coefficienti di amplificazione dei momenti flettenti
k_Minf.pl = 1; % fattore amplificativo del momento inferiore per avere armature maggiori
k_Msup.pl = 1;
k_Minf.tr = 1; % fattore amplificativo del momento inferiore per avere armature maggiori
k_Msup.tr = 1;

% passo di integrazione
dx = 0.001;

%% Dati armatura Platea

arm.pl.cf = 40;
arm.pl.template = [ arm.pl.cf, 0, 0;...
                    geom.pl.h - arm.pl.cf, 0, 0];  % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante

% armatura superiore
arm.pl.sup.fiv = 10:2:24;
arm.pl.sup.nb1 = 5;
arm.pl.sup.nb2 = 0;
arm.pl.sup.lock = 'true';

% armatura  inferiore
arm.pl.inf.fiv = 10:2:24;
arm.pl.inf.nb1 = 5;
arm.pl.inf.nb2 = 0;
arm.pl.inf.lock = 'false';

%% Dati armatura Trave.

arm.tr.cf = 50;
arm.tr.template = [ arm.tr.cf, 0, 0;...
                    280, 2, 14;...
                    520, 2, 14;...
                    geom.tr.h - arm.tr.cf,  0, 0]; % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante

% armatura superiore
arm.tr.sup.fiv = 14:2:26;
arm.tr.sup.nb1 = 4;
arm.tr.sup.nb2 = 1:4;
arm.tr.sup.lock = 'false';

% armatura inferiore
arm.tr.inf.fiv = 14:2:26;
arm.tr.inf.nb1 = 4;
arm.tr.inf.nb2 = 1:4;
arm.tr.inf.lock = 'false';

% dati armatura staffe per le travi
arm.tr.staffe.nb_sw = 4;
arm.tr.staffe.fi_sw = 8:2:14;
arm.tr.staffe.s_lim = [50, 10];
arm.tr.staffe.cf = 40; % copriferro


%% Funzione freccia elastica in mezzeria
f_el = @(q,l,E,J) -1/384 * q*l^4/(E*J);

%% Opzioni per il dimensionamento
opt1 = {'tipo', '''elastica''', 'precisione', '6'}; % per dimenSezione()
opt2 = {'piegate',true};    % per calcPesoCamp
