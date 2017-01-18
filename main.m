clearvars
clc
%% Dimensionamento Platea
% in condizioni SLE la platea è soggeta al carico q = 66.3 kN/m^2. Il
% criterio dimensionante è la freccia nel punto più sfavorevole, i.e. nel
% punto medio della platea compreso tra due pali. Pertanto la luce
% effettiva di calcolo è pari a (Lpl + Ltr)^2.

% freccia limite
f_lim = 1/2500;

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
% composizione del carico:
%   G1 = 16.3 kN/m2
%   Qk = 50.0 kN/m2
soll.pl.q = 50;  % [kN/m2] carico agente in condizioni SLE
soll.pl.qSLU = 96.19;  % [kN/m2] carico agente in condizioni SLU
soll.pl.N = 0;
soll.pl.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata
soll.pl.V = @(q,l,x) -q*x + q*l/2; % taglio lungo l'asse della trave
% dati sollecitazione della trave
soll.tr.q = 0;  % inizializzo il valore di carico da stimare per ogni iteraizione in funzione della luce di inflessione della platea
soll.tr.qSLU = 0;  % idem come sopra
soll.tr.N = 0;
soll.tr.M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave 1 campata
soll.tr.V = @(q,l,x) -q*x + q*l/2;

%% LIMITAZIONI AL PROBLEMA DI OTTIMIZZAZIONE
% Di seguito si elencano le limitazioni e le costanti che caratterizzano
% questo problema di ottimizzazione. Sono inserite le dimensioni della
% fondazione, i limiti per il numero di travi (i.e. luce libera di
% inflessione della platea, i limiti per il numero di pali/trave (i.e. luce
% libera di inflessione delle travi).
% In questa parte si indicano i costi unitari al peso ed al palo.
L.X = 56;   % [m] dimensione in x della fondazione
L.Y = 79;   % [m] dimensione in y della fondazione
L.b_max = 12; % [m] lunghezza massima delle barre
num.tr = 10:29; % numeri di travi nella fondazione 10:29
num.pali_tr = 14:40; % numero di pali per trave 14:40
cu.steel = 0.79; % costo unitatio acciaio [€/kg]
cu.palo  = 670.08;  % costo unitario palo [€/kg]
peso_acciaio = 7850;    % [kg/m^3] peso unitario dell'acciaio
RZ_palo.max = 960; % [kN] carico massimo ammissibile per il palo

%%
% costanti
cost.k = 4*sqrt(6)/9; % amplificatzione della freccia per ottenerre il massimo rapporto df/dL in funzione delle condizioni di appoggio (1 per app-app, 4*sqrt(6)/9 per inc-inc)
f_lim = f_lim/cost.k;
cost.c = 0.5;   % costante per il calcolo del parametro zeta (c = 0.5 per carichi ripetuti o lunga durata, c = 1 per carichi brevi)
cost.rot = 0;   % costante di rotazione (integrale della curvatura)
cost.spo = 0;   % costante di spostamento (integrale della rotazione)
% passo di integrazione
dx = 0.001;

%% Dati armatura Platea

arm.pl.template = [40, 0, 0; 260, 0, 0];  % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante

% armatura superiore
arm.pl.sup.fiv = 10:2:24;
arm.pl.sup.nb1 = 4:5;
arm.pl.sup.nb2 = 4:5;
arm.pl.sup.lock = 'true';

% armatura  inferiore
arm.pl.inf.fiv = 10:2:24;
arm.pl.inf.nb1 = 4:5;
arm.pl.inf.nb2 = 0;
arm.pl.inf.lock = 'false';

%% Dati armatura Trave.

arm.tr.template = [ 50, 0, 0;...
                    280, 2, 14;...
                    520, 2, 14;...
                    750,  0, 0]; % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante

% armatura superiore
arm.tr.sup.fiv = 14:2:24;
arm.tr.sup.nb1 = 4;
arm.tr.sup.nb2 = 1:4;
arm.tr.sup.lock = 'false';

% armatura inferiore
arm.tr.inf.fiv = 14:2:24;
arm.tr.inf.nb1 = 4;
arm.tr.inf.nb2 = 0;
arm.tr.inf.lock = 'false';

% dati armatura staffe per le travi
arm.tr.staffe.nb_sw = 4;
arm.tr.staffe.fi_sw = 10:2:14;
arm.tr.staffe.s_lim = [50, 10];
arm.tr.staffe.cf = 40; % copriferro


%% Funzione freccia elastica in mezzeria
f_el = @(q,l,E,J) -1/384 * q*l^4/(E*J);

%% Variazione luce libera di inflessione della platea
L.pl = L.X ./ (num.tr - 1);      % [m] variazione della luce libera di inflessione della platea
L.tr = L.Y ./ (num.pali_tr - 1); % [m] variazione della luce libera di inflessione della trave

%% Dati dei pali
num.pali = num.tr' * num.pali_tr;   % numero totale di pali nelle varie configurazioni size: lenght(num.tr) x length(numpali_tr)
RZ_palo.SLE = L.pl' * L.tr * soll.pl.q;
RZ_palo.SLU = L.pl' * L.tr * soll.pl.qSLU;
RZ_palo.cond = RZ_palo.SLE <= RZ_palo.max;

%% Inizio ciclo principale

% preallocamento variabili di carico
M.pl.inf = 0;
M.pl.sup = 0;
M.pl.mx = 0;
M.pl = repmat(M.pl, length(L.pl), 1);
M.tr.inf = 0;
M.tr.sup = 0;
M.tr.mx = 0;
M.tr = repmat(M.tr, length(L.pl), length(L.tr));
V.tr.mx = 0;
V.tr = repmat(V.tr, length(L.pl), length(L.tr));

% preallocamento risultati della freccia
ris.pl.s_min = ones(size(length(L.pl)))*nan;
ris.pl.s_rel = ones(size(length(L.pl)))*nan;
ris.pl.s_logico = ones(size(length(L.pl)))*nan;

ris.tr.s_min = ones(length(L.pl), length(L.tr))*nan;
ris.tr.s_rel = ones(length(L.pl), length(L.tr))*nan;
ris.tr.s_logico = ones(length(L.pl), length(L.tr))*nan;

ris.tot.s_min = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_rel = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_rel2 = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_logico = false(length(L.pl), length(L.tr));

soll.tr = repmat(soll.tr, length(L.pl), 1); % sollecitazioni per travi (funzione di L.pl)

% preallocamento struttura risultati dell'armatura
fieldname = {'nb1', 'nb2', 'fi1', 'fi2', 'As1', 'As2', 'As_tot','Mrd_fi1', 'Mrd', 'ratio', 'armatura', 'Lp_fi2', 'L_fi1', 'L_fi2', 'peso_fi1', 'peso_fi2', 'peso_tot'};
for i = 1:length(fieldname)
    if ~strcmp(fieldname{i},'Lp_fi2')
        ris.pl.arm_sup.(fieldname{i}) = nan;
    else
        ris.pl.arm_sup.(fieldname{i}) = struct;
    end
end
ris.pl.arm_sup = repmat(ris.pl.arm_sup, length(L.pl), 1);
ris.pl.arm_inf = ris.pl.arm_sup;

% è sufficiente ampliare ris.pl.arm_sup nell'altra direzione per ottenere
% una struttura NxM dove N = length(L.pl) e M = length(L.tr)
ris.tr.arm_sup = repmat(ris.pl.arm_sup, 1, length(L.tr));
ris.tr.arm_inf = ris.tr.arm_sup;

fieldname = {'nb', 'fi', 'Asw', 's', 's_max', 'Asw_s', 'cot_theta', 'Vrsd', 'Vrcd', 'Vrd', 'ratio', 'zona', 'L_zona', 'Q_st', 'Peso'};
for i = 1:length(fieldname)
    ris.tr.staffe.arm.(fieldname{i}) = nan;
end
ris.tr.staffe = repmat(ris.tr.staffe, length(L.pl), length(L.tr));

% preallocamento risultati del peso
ris.pl.peso.base = ones(length(L.pl),1) * nan;
ris.pl.peso.add = ones(length(L.pl),1) * nan;
ris.pl.peso.tot = ones(length(L.pl),1) * nan;

ris.tr.peso.base = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.add = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.staffe_tr = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.staffe = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.tot = ones(length(L.pl),length(L.tr)) * nan;

% preallocamento variabilli logiche
armaturaValida.pl = true;
armaturaValida.pl = repmat(armaturaValida.pl, length(L.pl), 1);
armaturaValida.tr = RZ_palo.cond;

% preallocamento risultati di costo
ris.costo.pl = ones(length(L.pl),1) * nan;
ris.costo.tr = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.pali = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.tot = ones(length(L.pl), length(L.tr)) * nan;

%% Calcolo della freccia
ite.current = 0;
ite.pl = 0;
ite.plMax = length(L.pl);
ite.trMax = length(L.tr);
ite.numMax = length(L.pl)*length(L.tr);
ite.progress = @(i) i/ite.numMax;
wb = waitbar(ite.progress(ite.current), sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));

%% corpo principale
for i_pl = 1:length(L.pl)
    
    %aggiornamento waitbar
    ite.pl = ite.pl + 1;
    ite.tr = 0;
    
    %% dimensionamento platea
    % Il processo di dimensionamento è composto dai seguenti step:
    %   1. Dimensionamento in semplice armatura per il lembo meno
    %   sollecitato;
    %   2. Ottimizzazione in funzione del peso dei risultati ottenuti;
    %   3. Aggiornamento dell'armatura tipo con i valori ottenuti (solo
    %   relativi alle barre 1
    %   4. Dimensionamento in armatura reale per la sezione maggiormente
    %   sollecitata
    %   5. Ottimizzazione in funzione del peso.
    
    M.pl(i_pl).inf = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
    M.pl(i_pl).sup = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0);    % momento all'appoggio, in funzione della luce
    M.pl(i_pl).mx = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0:dx:L.pl(i_pl));
    
    % dimensionamento dell'armatura superiore/inferiore
    
    if abs(M.pl(i_pl).inf) <= abs(M.pl(i_pl).sup)
        lembo = {'inf', 'sup'};
        l_sos = [length(arm.pl.template(:,1)), 1];
    else
        lembo = {'sup', 'inf'};
        l_sos = [1, length(arm.pl.template(:,1))];
    end
    
    %     for i_lem = 1:length(lembo)
    j = 0; % variabile di ciclo while
    continua_ciclo.sup = true;
    continua_ciclo.inf = true;
    while continua_ciclo.sup == true || continua_ciclo.inf == true
        i_lem = 1 + rem(j,2);
        % in funzione del lembo che si sta dimensionando, determino quale
        % lato deve essere aggiornato con i valori dell'iterazione
        % precedente
        
        % aggiornamento dell'armatura tipo
        if j == 0
            arm.pl.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [0, 0];    % per la prima iterazione, l'armatura è semplice
        else
            arm.pl.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [ris_.nb1, ris_.fi1];    % aggiorno con i valori dell'iterazione precedente
        end
        
        % calcolo momento resistente e peso
        ris_ = dimenSezione(geom.pl.sezione, arm.pl.template, arm.pl.(lembo{i_lem}).fiv, arm.pl.(lembo{i_lem}).nb1, arm.pl.(lembo{i_lem}).nb2, fck, M.pl(i_pl).(lembo{i_lem}), soll.pl.N, 'tipo', '''elastica''', 'lock', arm.pl.(lembo{i_lem}).lock, 'precisione', '6');
        ris_ = calcPesoCamp(ris_, sign(M.pl(i_pl).(lembo{i_lem})), L.pl(i_pl), M.pl(i_pl).mx, 'piegate', true);
        
        % questo ciclo while riduce le possibili soluzioni a quella con peso
        % minore
        while length(ris_) > 1
            if ris_(1).peso_tot < ris_(2).peso_tot
                ris_(2) = [];
            else
                ris_(1) = [];
            end
        end
        
        % CONDIZIONI DI FINE CICLO
        % se ris_ == [] significa che non c'è una configurazione di
        % armatura valida per quella sollecitazione, i.e. per L.pl
        if ~isempty(ris_)
            % condizione per terminare il ciclo: se nb1, nb2, fi1, fi2
            % calcolati nell'iterazione attuale (ris_) non sono cambiati
            % dall'iterazione precedente (ris), allora il ciclo può
            % interrompersi
            if j >= 2 && ...    % garantisce check delle condizioni a partire dalle 3a iterazione
                    ris_.nb1 == ris.pl.(['arm_' lembo{i_lem}])(i_pl).nb1 && ...
                    ris_.nb2 == ris.pl.(['arm_' lembo{i_lem}])(i_pl).nb2 && ...
                    ris_.fi1 == ris.pl.(['arm_' lembo{i_lem}])(i_pl).fi1 && ...
                    ris_.fi2 == ris.pl.(['arm_' lembo{i_lem}])(i_pl).fi2
                continua_ciclo.(lembo{i_lem}) = false;
            end
            ris.pl.(['arm_' lembo{i_lem}])(i_pl) = ris_; % salvo la soluzione temporanea nel struttura globale
        else
            armaturaValida.pl(i_pl) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
            break
        end
        j = j+1;
    end
    
    % variabili temporanee per il calcolo della freccia. Si considera solo
    % l'armatura di base per il calcolo della freccia perché in condizioni
    % semplici. La freccia è sovrastimata in favore di sicurezza;
    if armaturaValida.pl(i_pl)
        arm_.nb = [ris.pl.arm_inf(i_pl).armatura.nb1; ris.pl.arm_inf(i_pl).armatura.nb2];
        arm_.diam = [ris.pl.arm_inf(i_pl).armatura.fi1; ris.pl.arm_inf(i_pl).armatura.fi2];
        arm_.d = [ris.pl.arm_inf(i_pl).armatura.d; ris.pl.arm_inf(i_pl).armatura.d];
        % calcolo freccia
        ris.pl.s_min(i_pl) = calcoloFreccia(L.pl(i_pl), geom.pl, arm_, fck, cost.c, dx, soll.pl, 'semplificato', f_el);
    else
        ris.pl.s_min(i_pl) = nan;
    end
    
    %% dimensionamento della trave
    
    % aggiornamento variabile di carico agente sulla trave
    if armaturaValida.pl(i_pl)
        soll.tr(i_pl).q = soll.pl.q * L.pl(i_pl);
        soll.tr(i_pl).qSLU = soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
    else
        soll.tr(i_pl).q = nan;
        soll.tr(i_pl).qSLU = nan;
    end
    
    % loop principale per i_tr
    if armaturaValida.pl(i_pl)
        for i_tr = 1:length(L.tr)
            % debug
            if i_pl == 19 && i_tr == 4
                disp('yo')
            end
            
            %aggiornamento waitbar
            ite.tr = ite.tr + 1;
            ite.current = (ite.pl-1)*ite.trMax + ite.tr;
            waitbar(ite.progress(ite.current), wb, sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));
            set(wb, 'Name', sprintf('Progress \n%0.2f%%', ite.progress(ite.current)*100));
            
            % sollecitazione agente
            M.tr(i_pl, i_tr).inf = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), L.tr(i_tr)/2);    % momento in mezzeria, in funzione della luce
            M.tr(i_pl, i_tr).sup = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0);    % momento all'appoggio, in funzione della luce
            M.tr(i_pl, i_tr).mx = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0:dx:L.tr(i_tr));
            V.tr(i_pl, i_tr).vx = soll.tr(i_pl).V(soll.tr(i_pl).qSLU, L.tr(i_tr), 0:dx:L.tr(i_tr));
            
            % dimensionamento dell'armatura superiore/inferiore
            % dimensionamento dell'armatura superiore/inferiore
            
            if abs(M.tr(i_pl, i_tr).inf) <= abs(M.tr(i_pl, i_tr).sup)
                lembo = {'inf', 'sup'};
                l_sos = [length(arm.tr.template(:,1)), 1];
            else
                lembo = {'sup', 'inf'};
                l_sos = [1, length(arm.tr.template(:,1))];
            end
            
            %     for i_lem = 1:length(lembo)
            j = 0; % variabile di ciclo while
            continua_ciclo.sup = true;
            continua_ciclo.inf = true;
            while continua_ciclo.sup == true || continua_ciclo.inf == true
                i_lem = 1 + rem(j,2);
                
                % aggiornamento dell'armatura tipo
                if j == 0
                    arm.tr.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [0, 0];    % per la prima iterazione, l'armatura è semplice
                else
                    arm.tr.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [ris_.nb1, ris_.fi1];    % aggiorno con i valori dell'iterazione precedente
                end
                
                % calcolo momento resistente e peso
                ris_ = dimenSezione(geom.tr.sezione, arm.tr.template, arm.tr.(lembo{i_lem}).fiv, arm.tr.(lembo{i_lem}).nb1, arm.tr.(lembo{i_lem}).nb2, fck, M.tr(i_pl, i_tr).(lembo{i_lem}), soll.tr(i_pl).N, 'tipo', '''elastica''', 'lock', arm.tr.(lembo{i_lem}).lock, 'precisione', '6');
                ris_ = calcPesoCamp(ris_, sign(M.tr(i_pl, i_tr).(lembo{i_lem})), L.tr(i_tr), M.tr(i_pl, i_tr).mx, 'piegate', true);
                
                % questo ciclo while riduce le possibili soluzioni a quella con peso
                % minore
                while length(ris_) > 1
                    if ris_(1).peso_tot < ris_(2).peso_tot
                        ris_(2) = [];
                    else
                        ris_(1) = [];
                    end
                end
                % CONDIZIONI DI FINE CICLO
                % se ris_ == [] significa che non c'è una configurazione di
                % armatura valida per quella sollecitazione, i.e. per L.pl
                if ~isempty(ris_)
                    % condizione per terminare il ciclo: se nb1, nb2, fi1, fi2
                    % calcolati nell'iterazione attuale (ris_) non sono cambiati
                    % dall'iterazione precedente (ris), allora il ciclo può
                    % interrompersi
                    if j >= 2 && ...    % garantisce check delle condizioni a partire dalle 3a iterazione
                            ris_.nb1 == ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr).nb1 && ...
                            ris_.nb2 == ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr).nb2 && ...
                            ris_.fi1 == ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr).fi1 && ...
                            ris_.fi2 == ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr).fi2
                        continua_ciclo.(lembo{i_lem}) = false;
                    elseif j >= 4 && ... % garantisce check delle condizioni a partire dalle 5a iterazione
                            ris_.nb1 == ris2_.(['arm_' lembo{i_lem}]).nb1 && ...
                            ris_.nb2 == ris2_.(['arm_' lembo{i_lem}]).nb2 && ...
                            ris_.fi1 == ris2_.(['arm_' lembo{i_lem}]).fi1 && ...
                            ris_.fi2 == ris2_.(['arm_' lembo{i_lem}]).fi2
                        continua_ciclo.(lembo{i_lem}) = false;
                    else
                        continua_ciclo.(lembo{i_lem}) = true;
                    end
                    ris2_.(['arm_' lembo{i_lem}]) = ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr); % valore dell'iterazione i-2
                    ris.tr.(['arm_' lembo{i_lem}])(i_pl, i_tr) = ris_; % salvo la soluzione temporanea nel struttura globale
                else
                    armaturaValida.tr(i_pl, i_tr) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
                    break
                end
                j = j+1;
            end
            
            % calcolo dell'armatura a taglio
            if armaturaValida.tr(i_pl, i_tr)
                [armSt_, pesoSt_]  = calcPesoStaffe(geom.tr.sezione, L.tr(i_tr)*1e3, arm.tr.template(:,1), arm.tr.staffe.nb_sw, arm.tr.staffe.fi_sw, ris_.fi1, arm.tr.staffe.s_lim, fck, soll.tr(i_pl).N, V.tr(i_pl, i_tr).vx, 'copriferro', arm.tr.staffe.cf);
                if ~isempty(armSt_)
                    ris.tr.staffe(i_pl, i_tr).arm = armSt_;
                    ris.tr.peso.staffe_tr(i_pl, i_tr) = pesoSt_;    % peso delle staffe per concio di trave
                else
                    armaturaValida.tr(i_pl, i_tr) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
                end
            end
            
            % variabili temporanee per il calcolo della freccia. Si
            % considera solo l'armatura di base per il calcolo della
            % freccia perché siamo in condizioni semplici. La freccia è
            % sovrastimata in favore di sicurezza;
            if armaturaValida.tr(i_pl, i_tr)
                arm_.nb = [ris.tr.arm_inf(i_pl, i_tr).armatura.nb1; ris.tr.arm_inf(i_pl, i_tr).armatura.nb2];
                arm_.diam = [ris.tr.arm_inf(i_pl, i_tr).armatura.fi1; ris.tr.arm_inf(i_pl, i_tr).armatura.fi2];
                arm_.d = [ris.tr.arm_inf(i_pl, i_tr).armatura.d; ris.tr.arm_inf(i_pl, i_tr).armatura.d];
                % calcolo freccia
                ris.tr.s_min(i_pl, i_tr) = calcoloFreccia(L.tr(i_tr), geom.tr, arm_, fck, cost.c, dx, soll.tr(i_pl), 'semplificato', f_el);
            else
                ris.tr.s_min(i_pl, i_tr) = nan;
            end
        end
    end
end
close(wb)
clear('h')

%% Combinazione dei risultati
for i_pl = 1:length(L.pl)
    % risultati per la platea
    ris.pl.s_rel(i_pl) = abs(ris.pl.s_min(i_pl))/(L.pl(i_pl)*1e3);
    ris.pl.s_logico(i_pl) =  ris.pl.s_rel(i_pl) <= f_lim;
    
    for i_tr = 1:length(L.tr)
        if and(armaturaValida.pl(i_pl), armaturaValida.tr(i_pl, i_tr))
            % risultati per le travi
            ris.tr.s_rel(i_pl, i_tr) = abs(ris.tr.s_min(i_pl, i_tr)/(L.tr(i_tr)*1e3));
            ris.tr.s_logico(i_pl, i_tr) = ris.tr.s_rel(i_pl, i_tr) <= f_lim;
            
            % platea + travi
            L.tot(i_pl, i_tr) = sqrt(L.pl(i_pl)^2 + L.tr(i_tr)^2);  % lunghezza combinata su cui verificare la freccia limite
            ris.tot.s_min(i_pl, i_tr) = ris.pl.s_min(i_pl) + ris.tr.s_min(i_pl, i_tr);
            ris.tot.s_rel(i_pl, i_tr) = abs(ris.tot.s_min(i_pl, i_tr)/(L.tot(i_pl, i_tr)*1e3));    % freccia in termini di luce libera di inflessione
            ris.tot.s_rel2(i_pl, i_tr) = abs(ris.pl.s_min(i_pl)*1e-3/L.pl(i_pl) + ris.tr.s_min(i_pl, i_tr)*1e-3/L.tr(i_tr));
            ris.tot.s_logico(i_pl, i_tr) = ris.tot.s_rel(i_pl, i_tr) <= f_lim;
        end
    end
end

%% Ottimizzazione dei risultati accettabili in funzione del costo
for i_pl = 1:length(L.pl)
    if sum(ris.tot.s_logico(i_pl, :))   % almeno un componenente è pari ad 1, pertanto vale la pena eseguire il calcolo
        clear('l_.pl')
        % calcolo del peso delle barre della platea
        l_.pl = struct;
        % Combina l'armatura superiore e quella inferiore in una struttura
        % unica. Usa l'armatura inferiore come template e sostituisce la
        % prima riga con il valori corrispondente dell'armatura superiore e
        % aggiorna la variabile di peso.
        ris.pl.comb_arm(i_pl) = ris.pl.arm_inf(i_pl).armatura;
        fieldname_ = fields( ris.pl.arm_sup(i_pl).armatura );
        for i = 1:length( fieldname_)-1 % escludo il campo 'Peso' che è l'ultimo
            ris.pl.comb_arm(i_pl).(fieldname_{i})(1) = ris.pl.arm_sup(i_pl).armatura.(fieldname_{i})(1);
        end
        ris.pl.comb_arm(i_pl).Peso = sum( ris.pl.comb_arm(i_pl).peso_tot ); % aggiorno il campo 'Peso'
        
        % determino qual è il massimo diametro per l'armatura continua
        % (nb1, fi1) e calcolo la lunghezza di ancoraggio.
        l_.pl.anc = 60 * ris.pl.comb_arm(i_pl).fi1*1e-3;    % [m] vettore delle lunghezze di ancoraggio in scarsa aderenza per barre continue
        l_.pl.tra = 0.9 * max(ris.pl.comb_arm(i_pl).d)/2*1e-3;   % [m] valore della lunghezza di traslazione dei momenti (uguale per tuttle barre)
        l_.pl.beff = L.b_max - l_.pl.anc - l_.pl.tra;   % [m] lunghezza effettiva delle barre
        l_.pl.k_fi1 = L.b_max ./ l_.pl.beff;   % fattore amplificativo che tiene conto delle sovrapposizione delle barre
        
        % calcola il peso delle barre fi1 i.e. l'armatura di base continua.
        % Il campo "peso_fi1" è moltiplicato per il coefficiente k_fi1 e
        % per il numero di campate della platea, pari al numero di travi -
        % 1. Questa quantità deve quindi essere moltiplicata per la
        % dimensione L.Y/b, dove b è la larghezza della sezione di calcolo,
        % i.e. n = 1 m;
        ris.pl.peso.base(i_pl) = sum(l_.pl.k_fi1 .* ris.pl.comb_arm(i_pl).peso_fi1) * ( num.tr(i_pl)-1 ) * L.Y/(geom.pl.b*1e-3);
        
        % calcola il peso delle barre addizionali, non continue,
        % pertanto il peso è pari al peso_fi2.
        ris.pl.peso.add(i_pl) = sum(ris.pl.comb_arm(i_pl).peso_fi2) * (num.tr(i_pl)-1) * L.Y/(geom.pl.b * 1e-3);
        
        % Peso totale della platea in funzione di L.pl
        ris.pl.peso.tot(i_pl) = ris.pl.peso.base(i_pl) + ris.pl.peso.add(i_pl);
        
        % calcolo del peso delle barre delle travi
        l_.tr = struct;
        for i_tr = 1:length(L.tr)
            if ris.tot.s_logico(i_pl, i_tr)
                clear('l_.tr')
                % calcolo del peso delle barre della trave
                l_.tr = struct;
                
                % Combina l'armatura superiore e quella inferiore in una struttura
                % unica. Usa l'armatura inferiore come template e sostituisce la
                % prima riga con il valori corrispondente dell'armatura superiore e
                % aggiorna la variabile di peso.
                ris.tr.comb_arm(i_pl, i_tr) = ris.tr.arm_inf(i_pl, i_tr).armatura;
                fieldname_ = fields( ris.tr.arm_sup(i_pl, i_tr).armatura );
                for i = 1:length( fieldname_)-1 % escludo il campo 'Peso' che è l'ultimo
                    ris.tr.comb_arm(i_pl, i_tr).(fieldname_{i})(1) = ris.tr.arm_sup(i_pl, i_tr).armatura.(fieldname_{i})(1);
                end
                ris.tr.comb_arm(i_pl, i_tr).Peso = sum( ris.tr.comb_arm(i_pl, i_tr).peso_tot ); % aggiorno il campo 'Peso'
                
                % determino qual è il massimo diametro per l'armatura continua
                % (nb1, fi1) e calcolo la lunghezza di ancoraggio.
                l_.tr.anc = 60 * ris.tr.comb_arm(i_pl, i_tr).fi1*1e-3;    % [m] vettore delle lunghezze di ancoraggio in scarsa aderenza per barre continue
                l_.tr.tra = 0.9 * max(ris.tr.comb_arm(i_pl, i_tr).d)/2*1e-3;   % [m] valore della lunghezza di traslazione dei momenti (uguale per tuttle barre)
                l_.tr.beff = L.b_max - l_.tr.anc - l_.tr.tra;   % [m] lunghezza effettiva delle barre
                l_.tr.k_fi1 = L.b_max ./ l_.tr.beff;   % fattore amplificativo che tiene conto delle sovrapposizione delle barre
                
                % calcola il peso delle barre fi1 i.e. l'armatura di base continua.
                % Il campo "peso_fi1" è moltiplicato per il coefficiente k_fi1 e
                % per il numero di campate della platea, pari al numero di travi -
                % 1. Questa quantità deve quindi essere moltiplicata per la
                % dimensione L.Y/b, dove b è la larghezza della sezione di calcolo,
                % i.e. n = 1 m;
                ris.tr.peso.base(i_pl, i_tr) = sum(l_.tr.k_fi1 .* ris.tr.comb_arm(i_pl, i_tr).peso_fi1) * ( num.pali_tr(i_tr)-1 ) * num.tr(i_pl);
                
                % calcola il peso delle barre addizionali, non continue,
                % pertanto il peso è pari al peso_fi2.
                ris.tr.peso.add(i_pl, i_tr) = sum(ris.tr.comb_arm(i_pl, i_tr).peso_fi2) * (num.pali_tr(i_tr)-1) * num.tr(i_pl);
                
                % calcolo del peso delle staffe
                ris.tr.peso.staffe(i_pl, i_tr) = ris.tr.peso.staffe_tr(i_pl, i_tr) * num.tr(i_pl) * (num.pali_tr(i_tr) - 1);   
                
                % Peso totale della trave in funzione di L.tr
                ris.tr.peso.tot(i_pl, i_tr) = ris.tr.peso.base(i_pl, i_tr) + ris.tr.peso.add(i_pl, i_tr) + ris.tr.peso.staffe(i_pl, i_tr);
            end
        end
    end
end
%% Computo del costo
for i_pl = 1:length(L.pl)
    if sum(ris.tot.s_logico(i_pl, :))
        % costo delle platea in funzione di L.pl
        ris.costo.pl(i_pl) = ris.pl.peso.tot(i_pl) * cu.steel;
        for i_tr = 1:length(L.tr)
            if ris.tot.s_logico(i_pl, i_tr)
                ris.costo.tr(i_pl, i_tr) = ris.tr.peso.tot(i_pl, i_tr) * cu.steel;
                ris.costo.pali(i_pl, i_tr) = num.pali(i_pl, i_tr) * cu.palo;
                ris.costo.tot(i_pl, i_tr) = ris.costo.pl(i_pl) + ris.costo.tr(i_pl, i_tr) + ris.costo.pali(i_pl, i_tr);
            end
        end
    end
end

[costo_min, I] = min(ris.costo.tot(:));
[i_min.pl, i_min.tr] = ind2sub(size(ris.costo.tot), I);
%% grafico delle frecce
figure(1)
subplot(1,2,1)
surf(L.pl, L.tr, ris.tot.s_rel')
title('Deformazioni massime')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('(f_{pl}+f_{tr})/L_{tot}')
hold on
surf(L.pl, L.tr, f_lim*ones(length(L.pl), length(L.tr))')
subplot(1,2,2)
contourf(L.pl, L.tr, ris.tot.s_logico',[0, 1]);
title('Soluzioni Valide')
xlabel('L_{platea}')
ylabel('L_{trave}')
%% grafico dei costi
figure(2)

subplot(3,2,1)
plot(L.pl, ris.costo.pl)
title('Costo Platea')
xlabel('L_{platea}')
ylabel('Costo Platea')

subplot(3,2,3)
[C.tr, h.tr] = contourf(L.pl, L.tr, ris.costo.tr');
clabel(C.tr, h.tr)
title('Costo Travi')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Costo Travi')

subplot(3,2,5)
[C.pali, h.pali] = contourf(L.pl, L.tr, ris.costo.pali');
clabel(C.pali, h.pali);
title('Costo Pali')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Costo pali')

subplot(3,2,[2, 4, 6])
[C.tot, h.tot] = contourf(L.pl, L.tr, ris.costo.tot');
clabel(C.tot, h.tot);
title('Costo Totale [%]')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Somma Costi')

%%
clc
out = fopen('results.txt','w+');
fprintf(out, 'La soluzione ottimale si ha in corripondenza degli indici:\n');
fprintf(out, '\t% 4s\t% 8s\t% 8s\n', '', 'min', 'L[m]');
label = {'pl', 'tr'};
for i = 1:length(label)
    fprintf(out, '\ti_%2s\t% 8.0f\t% 8.2f \n', label{i}, i_min.(label{i}), L.(label{i})(i_min.(label{i})));
end
fprintf(out, '\nChe corrispondono a %d travi e %d pali per trave, per un totale di %d pali.\n', num.tr(i_min.pl), num.pali_tr(i_min.tr), num.pali(i_min.pl, i_min.tr));
fprintf(out, '\nLe deformazioni minime per questa combinazione sono pari a:\n');
fprintf(out, '\t% -4s\t% 12s\t% 12s\t% s<%.2e\n', '', 's_min[mm]', 's_rel[%]', 's_rel', f_lim );
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'pl', ris.pl.s_min(i_min.pl), ris.pl.s_rel(i_min.pl), ris.pl.s_logico(i_min.pl));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tr', ris.tr.s_min(i_min.pl, i_min.tr), ris.tr.s_rel(i_min.pl, i_min.tr), ris.tr.s_logico(i_min.pl, i_min.tr));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot1', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel(i_min.pl, i_min.tr), ris.tot.s_logico(i_min.pl, i_min.tr));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot2', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr)<= f_lim);
fprintf(out, '\nIl costo totale stimato è pari a %.2f €', ris.costo.tot(i_min.pl, i_min.tr));
fclose(out);
type('results.txt')