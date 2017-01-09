clearvars
clc
%% Dimensionamento Platea
% in condizioni SLE la platea è soggeta al carico q = 66.3 kN/m^2. Il
% criterio dimensionante è la freccia nel punto più sfavorevole, i.e. nel
% punto medio della platea compreso tra due pali. Pertanto la luce
% effettiva di calcolo è pari a (Lpl + Ltr)^2.

% freccia limite
f_lim = 1/3000;

% Dati geometrici sezione della platea
geom.pl.b = 1000;
geom.pl.h = 300;
geom.pl.x0 = 0;
geom.pl.y0 = 0;
geom.pl.sezione = rettangolo(geom.pl.b, geom.pl.h, geom.pl.x0, geom.pl.y0, 1, 100); % preparazione della sezione discretizzata

% Dati geometrici sezione della trave
geom.tr.b = 500;
geom.tr.h = 800;
geom.tr.x0 = 0;
geom.tr.y0 = 0;
geom.tr.sezione = rettangolo(geom.tr.b, geom.tr.h, geom.tr.x0, geom.tr.y0, 1, 100); % preparazione della sezione discretizzata


% dati materiale
fck = 25;

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
cu.steel = 0.80; % costo unitatio acciaio [€/kg]
cu.palo  = 600;  % costo unitario palo [€/kg]
peso_acciaio = 7850;    % [kg/m^3] peso unitario dell'acciaio

%%
% costanti
cost.c = 0.5;
cost.rot = 0;
cost.spo = 0;
supinf = {'sup', 'inf'};    % per calcolo armatura superiore/inferiore
% passo di integrazione
dx = 0.001;
% dati armatura platea tipo
% arm superiore
arm.pl.sup.fi_lim = [10, 24];
arm.pl.sup.nb1 = 5;
arm.pl.sup.nb2 = 5;
arm.pl.sup.lock = 'yes';
arm.pl.sup.d = 40;
% arm inferiore
arm.pl.inf.fi_lim = [10, 24];
arm.pl.inf.nb1 = 5;
arm.pl.inf.nb2 = 0;
arm.pl.inf.lock = 'no';
arm.pl.inf.d = 260;

% dati armatura trave tipo
arm.tr.fi_lim = [10, 24];
arm.tr.nb1 = [3, 5];
arm.tr.nb2 = [0, 4];
arm.tr.lock = 'no';
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
%arm.pl = repmat(arm.pl, length(L.pl),1);   % prealloco la struttura arm.pl (N)
%arm.tr = repmat(arm.tr, length(L.pl), length(L.tr)); % prealloco la struttura arm.tr  (NxM)
M.pl.max = 0;
M.pl.min = 0;
M.pl.mx = 0;
M.pl = repmat(M.pl, length(L.pl), 1);
M.tr.max = 0;
M.tr.min = 0;
M.tr.mx = 0;
M.tr = repmat(M.tr, length(L.pl), length(L.tr));

ris.pl.s_min = zeros(size(length(L.pl)));
ris.tr.s_min = zeros(length(L.pl), length(L.tr));
soll.tr = repmat(soll.tr, length(L.pl), 1);
%%
for i_pl = 1:length(L.pl)
    % dimensionamento platea
    M.pl(i_pl).inf = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
    M.pl(i_pl).sup = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0);    % momento all'appoggio, in funzione della luce
    M.pl(i_pl).mx = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0:dx:L.pl);
    
    % dimensionamento dell'armatura superiore/inferiore
    armaturaValida.pl(i_pl) = true;
    for i_si = 1:length(supinf)
        sezione_ = geom.pl.sezione;
        d_ = arm.pl.(supinf{i_si}).d;
        % se si calcola il momento resistente negativo è necessario
        % invertire i parametri inerenti l'altezza della sezione
        % (inversione dell'asse delle y)
        if strcmp(supinf{i_si},'sup')
            sezione_(:,2) = geom.pl.h - geom.pl.sezione(:,2);
            d_ = geom.pl.h - arm.pl.(supinf{i_si}).d;
        end
        armTemp = dimenSezione(sezione_, d_, arm.pl.(supinf{i_si}).fi_lim, arm.pl.(supinf{i_si}).nb1, arm.pl.(supinf{i_si}).nb2, fck, M.pl(i_pl).(supinf{i_si}), soll.pl.N, 'tipo', 'elastica', 'lock', arm.pl.(supinf{i_si}).lock, 'precisione', 6);
        armTemp = calcPesoCamp(armTemp, d_, sign(M.pl(i_pl).(supinf{i_si})), L.pl(i_pl), M.pl(i_pl).mx, 'piegate', 'yes');
        % questo ciclo while riduce le possibili soluzioni a quella con peso
        % minore
        while length(armTemp) > 1
            if armTemp(1).peso_tot < armTemp(2).peso_tot
                armTemp(2) = [];
            else
                armTemp(1) = [];
            end
        end
        % se armTemp == [] significa che non c'è una configurazione di
        % armatura valida per quella sollecitazione, i.e. per L.pl
        if ~isempty(armTemp)   
            ris.pl.arm.(supinf{i_si})(i_pl) = armTemp; % salvo la soluzione temporanea nel struttura globale
        else
            armaturaValida.pl(i_pl) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
        end
    end
    % variabili temporanee per il calcolo della freccia. Si considera solo
    % l'armatura di base per il calcolo della freccia perché in condizioni
    % semplici. La freccia è sovrastimata in favore di sicurezza;
    if armaturaValida.pl(i_pl)
        arm_.nb = [ris.pl.arm.sup(i_pl).nb1; ris.pl.arm.inf(i_pl).nb1];
        arm_.diam = [ris.pl.arm.sup(i_pl).fi1; ris.pl.arm.inf(i_pl).fi1];
        arm_.d = [arm.pl.sup.d; arm.pl.inf.d];
        % calcolo freccia
        ris.pl.s_min(i_pl) = calcoloFreccia(L.pl(i_pl), geom.pl, arm_, fck, cost.c, soll.pl, dx, 'semplificato', f_el);
    else
        ris.pl.s_min(i_pl) = nan;
    end
            
%     % aggiornamento variabile di carico agente sulla trave
%     soll.tr(i_pl).q = soll.pl.q * L.pl(i_pl);
%     soll.tr(i_pl).qSLU = soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
%     
%     % dimensionamento della trave
%     for i_tr = 1:length(L.tr)        
%         M_max.tr(i_pl, i_tr) = soll.tr(i_pl).M( soll.tr(i_pl).qSLU, L.tr(i_tr), L.tr(i_tr)/2);
%         As_min = M_max.tr(i_pl, i_tr)*1e6 / (0.9*max(arm.tr(i_pl, i_tr).d) * mat.steel.fyd) /arm.tr(i_pl, i_tr).nb(2);   % armatura minima con criterio semplificato dimensionata allo SLU
%         diam_min = ceil( ceil( sqrt( 4*As_min/pi ) ) /2 ) * 2;   % arrotonda ad un diametro reale
%         arm.tr(i_pl, i_tr).diam = [diam_min;14;14;diam_min];
%         % calcolo freccia della trave
%         ris.tr.s_min(i_pl, i_tr) = calcoloFreccia(L.tr(i_tr), geom.tr, arm.tr(i_pl, i_tr), fck, cost.c, dx, soll.tr(i_pl), 'semplificato', f_el);
%         
%     end
end
% %% Combinazione dei risultati
% for i_pl = 1:length(L.pl)
%     for i_tr = 1:length(L.tr)
%         L.tot(i_pl, i_tr) = sqrt(L.pl(i_pl)^2 + L.tr(i_tr)^2);  % lunghezza combinata su cui verificare la freccia limite
%         ris.tot.s_min(i_pl, i_tr) = ris.pl.s_min(i_pl) + ris.tr.s_min(i_pl, i_tr);
%         ris.tot.s_rel(i_pl, i_tr) = ris.tot.s_min(i_pl, i_tr)/(L.tot(i_pl, i_tr)*1e3);    % freccia in termini di luce libera di inflessione      
%     end
% end
