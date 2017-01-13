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
num.tr = 10:30; % numeri di travi nella fondazione
num.pali_tr = 14:40; % numero di pali per trave
cu.steel = 0.79; % costo unitatio acciaio [€/kg]
cu.palo  = 670.08;  % costo unitario palo [€/kg]
peso_acciaio = 7850;    % [kg/m^3] peso unitario dell'acciaio

%%
% costanti
cost.c = 0.5;
cost.rot = 0;
cost.spo = 0;
supinf = {'sup', 'inf'};    % per calcolo armatura superiore/inferiore
% passo di integrazione
dx = 0.001;

% dati armatura platea superiore
arm.pl.sup.fi_lim = [10, 24];
arm.pl.sup.nb1 = 3:5;
arm.pl.sup.nb2 = 5;
arm.pl.sup.lock = 'yes';
arm.pl.sup.d = 40;

% arm inferiore platea inferiore
arm.pl.inf.fi_lim = [10, 24];
arm.pl.inf.nb1 = 3:5;
arm.pl.inf.nb2 = 0;
arm.pl.inf.lock = 'no';
arm.pl.inf.d = 260;

% dati armatura trave superiore
arm.tr.sup.fi_lim = [14, 24];
arm.tr.sup.nb1 = [3, 5];
arm.tr.sup.nb2 = [0, 4];
arm.tr.sup.lock = 'no';
arm.tr.sup.d = 50;

% dati armatura trave inferiore
% barre intermedie di costruzione
arm.tr.med.nb = [2; 2];
arm.tr.med.fi = [14; 14];
arm.tr.med.d = [270; 580];

% dati armatura trave inferiore
arm.tr.inf.fi_lim = [14, 24];
arm.tr.inf.nb1 = [3, 5];
arm.tr.inf.nb2 = [0, 4];
arm.tr.inf.lock = 'no';
arm.tr.inf.d = 750;

%% Funzione freccia elastica in mezzeria
f_el = @(q,l,E,J) -1/384 * q*l^4/(E*J);

%% Variazione luce libera di inflessione della platea
L.pl = L.X ./ (num.tr - 1);      % [m] variazione della luce libera di inflessione della platea
L.tr = L.Y ./ (num.pali_tr - 1); % [m] variazione della luce libera di inflessione della trave

%%
num.pali = num.tr' * num.pali_tr;   % numero totale di pali nelle varie configurazioni size: lenght(num.tr) x length(numpali_tr)



%% Inizio ciclo principale

% preallocamento variabili di carico
M.pl.max = 0;
M.pl.min = 0;
M.pl.mx = 0;
M.pl = repmat(M.pl, length(L.pl), 1);
M.tr.max = 0;
M.tr.min = 0;
M.tr.mx = 0;
M.tr = repmat(M.tr, length(L.pl), length(L.tr));

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
fieldname = {'nb1', 'nb2', 'fi1', 'fi2', 'As1', 'As2', 'As_tot','Mrd_fi1', 'Mrd', 'ratio', 'Lp_fi2', 'L_fi1', 'L_fi2', 'peso_fi1', 'peso_fi2', 'peso_tot'};
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

% preallocamento risultati del peso
ris.pl.peso.base = ones(length(L.pl),1) * nan;
ris.pl.peso.add = ones(length(L.pl),1) * nan;
ris.pl.peso.tot = ones(length(L.pl),1) * nan;

ris.tr.peso.base = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.add = ones(length(L.pl),length(L.tr)) * nan;
ris.tr.peso.tot = ones(length(L.pl),length(L.tr)) * nan;

% preallocamento variabilli logiche
armaturaValida.pl = true;
armaturaValida.pl = repmat(armaturaValida.pl, length(L.pl), 1);
armaturaValida.tr = repmat(armaturaValida.pl, 1, length(L.tr));

% preallocamento risultati di costo
ris.costo.pl = ones(length(L.pl),1) * nan;
ris.costo.tr = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.pali = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.tot = ones(length(L.pl), length(L.tr)) * nan;

%% Calcolo della freccia
ite.current = 0;
ite.numMax = length(L.pl)*length(L.tr);
ite.progress = @(i) i/ite.numMax;
hw = waitbar(ite.progress(ite.current), sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));
%%
for i_pl = 1:length(L.pl)
    
    %% dimensionamento platea
    M.pl(i_pl).inf = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
    M.pl(i_pl).sup = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0);    % momento all'appoggio, in funzione della luce
    M.pl(i_pl).mx = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0:dx:L.pl(i_pl));
    
    % dimensionamento dell'armatura superiore/inferiore
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
            ris.pl.(['arm_' supinf{i_si}])(i_pl) = armTemp; % salvo la soluzione temporanea nel struttura globale
        else
            armaturaValida.pl(i_pl) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
            break
        end
    end
    
    % variabili temporanee per il calcolo della freccia. Si considera solo
    % l'armatura di base per il calcolo della freccia perché in condizioni
    % semplici. La freccia è sovrastimata in favore di sicurezza;
    if armaturaValida.pl(i_pl)
        arm_.nb = [ris.pl.arm_sup(i_pl).nb1; ris.pl.arm_inf(i_pl).nb1];
        arm_.diam = [ris.pl.arm_sup(i_pl).fi1; ris.pl.arm_inf(i_pl).fi1];
        arm_.d = [arm.pl.sup.d; arm.pl.inf.d];
        % calcolo freccia
        ris.pl.s_min(i_pl) = calcoloFreccia(L.pl(i_pl), geom.pl, arm_, fck, cost.c, dx, soll.pl, 'semplificato', f_el);
    else
        ris.pl.s_min(i_pl) = nan;
    end
    
    % aggiornamento variabile di carico agente sulla trave
    if armaturaValida.pl(i_pl)
        soll.tr(i_pl).q = soll.pl.q * L.pl(i_pl);
        soll.tr(i_pl).qSLU = soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
    else
        soll.tr(i_pl).q = nan;
        soll.tr(i_pl).qSLU = nan;
    end
    
    
    %% dimensionamento della trave
    for i_tr = 1:length(L.tr)
        ite.current = ite.current +1;
        waitbar(ite.progress(ite.current), hw, sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));
        if armaturaValida.pl(i_pl)
            % sollecitazione agente
            M.tr(i_pl, i_tr).inf = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), L.tr(i_tr)/2);    % momento in mezzeria, in funzione della luce
            M.tr(i_pl, i_tr).sup = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0);    % momento all'appoggio, in funzione della luce
            M.tr(i_pl, i_tr).mx = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0:dx:L.tr(i_tr));
            
            % dimensionamento dell'armatura superiore/inferiore
            for i_si = 1:length(supinf)
                sezione_ = geom.tr.sezione;
                d_ = arm.tr.(supinf{i_si}).d;
                % se si calcola il momento resistente negativo è necessario
                % invertire i parametri inerenti l'altezza della sezione
                % (inversione dell'asse delle y)
                if strcmp(supinf{i_si},'sup')
                    sezione_(:,2) = geom.tr.h - geom.tr.sezione(:,2);
                    d_ = geom.tr.h - arm.tr.(supinf{i_si}).d;
                end
                armTemp = dimenSezione(sezione_, d_, arm.tr.(supinf{i_si}).fi_lim, arm.tr.(supinf{i_si}).nb1, arm.tr.(supinf{i_si}).nb2, fck, M.tr(i_pl, i_tr).(supinf{i_si}), soll.tr(i_pl).N, 'tipo', 'elastica', 'lock', arm.tr.(supinf{i_si}).lock, 'precisione', 6);
                armTemp = calcPesoCamp(armTemp, d_, sign(M.tr(i_pl, i_tr).(supinf{i_si})), L.tr(i_tr), M.tr(i_pl, i_tr).mx, 'piegate', 'yes');
                
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
                    ris.tr.(['arm_' supinf{i_si}])(i_pl, i_tr) = armTemp; % salvo la soluzione temporanea nel struttura globale
                else
                    armaturaValida.tr(i_pl, i_tr) = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
                    break
                end
            end
            
            % variabili temporanee per il calcolo della freccia. Si considera solo
            % l'armatura di base per il calcolo della freccia perché in condizioni
            % semplici. La freccia è sovrastimata in favore di sicurezza;
            if armaturaValida.tr(i_pl, i_tr)
                arm_.nb = [ris.tr.arm_sup(i_pl, i_tr).nb1; arm.tr.med.nb; ris.tr.arm_inf(i_pl, i_tr).nb1];
                arm_.diam = [ris.tr.arm_sup(i_pl, i_tr).fi1; arm.tr.med.fi; ris.tr.arm_inf(i_pl, i_tr).fi1];
                arm_.d = [arm.tr.sup.d; arm.tr.med.d; arm.tr.inf.d];
                % calcolo freccia
                ris.tr.s_min(i_pl, i_tr) = calcoloFreccia(L.tr(i_tr), geom.tr, arm_, fck, cost.c, dx, soll.tr(i_pl), 'semplificato', f_el);
            else
                ris.tr.s_min(i_pl, i_tr) = nan;
            end
        end
    end
end
close(hw)
clear('hw')
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
        ris.pl.peso.base(i_pl) = 0;
        ris.pl.peso.add(i_pl) = 0;
        l_.pl = struct;
        for i_si = 1:length(supinf)
            % i seguenti parametri sono necessari per il calcolo
            % dell'armatura continua su tutta la lunghezza di L.X
            l_.pl.anc.(supinf{i_si}) = 60 * ris.pl.(['arm_' supinf{i_si}])(i_pl).fi1 * 1e-3;   % [m]lunghezza di ancoragggio per il diametro fi1 nell'ipotesi scarsa aderenza
            l_.pl.tra.(supinf{i_si})= 0.9 * ((2-i_si)*geom.pl.h + (2*(i_si==2) - 1) * arm.pl.(supinf{i_si}).d) * 1e-3; % [m] lunghezza di traslazione del momento
            l_.pl.beff.(supinf{i_si}) = L.b_max - l_.pl.anc.(supinf{i_si}) - l_.pl.tra.(supinf{i_si});   % [m] lunghezza effettiva delle barre
            l_.pl.Xeff.(supinf{i_si}) = L.X/l_.pl.beff.(supinf{i_si}) * L.b_max;   % [m] lunghezza totale delle barre > di L.X perché considera anche le sovrapposizioni
            % calcola il peso delle barre fi1 i.e. l'armatura di base
            % continua. Il campo "peso_fi1" è esrpresso in [kg] pertanto è
            % necessario dividerlo per la lunghezza della campata per
            % ottenere l'incidenza al metro lineare che moltiplica la
            % lunghezza efficace l_.pl.Xeff.(...).
            ris.pl.peso.base(i_pl) = ris.pl.peso.base(i_pl) + l_.pl.Xeff.(supinf{i_si}) * ris.pl.(['arm_' supinf{i_si}])(i_pl).peso_fi1/ris.pl.(['arm_' supinf{i_si}])(i_pl).L_fi1 * L.Y/(geom.pl.b*1e-3);   % il peso è moltiplicato per L.Y/geom.pl.b per ottenere il totale esteso a tutta la superficie della platea
            % calcola il peso delle barre addizionali, non continue,
            % pertanto il peso è pari al peso_fi2
            ris.pl.peso.add(i_pl) = ris.pl.peso.add(i_pl) + ris.pl.(['arm_' supinf{i_si}])(i_pl).peso_fi2 * (num.tr(i_pl)-1) * L.Y/(geom.pl.b * 1e-3);
            % Peso totale della platea in funzione di L.pl
            ris.pl.peso.tot(i_pl) = ris.pl.peso.base(i_pl) + ris.pl.peso.add(i_pl);
        end
        
        % calcolo del peso delle barre delle travi
        l_.tr = struct;
        for i_tr = 1:length(L.tr)
            ris.tr.peso.base(i_pl, i_tr) = 0;
            ris.tr.peso.add(i_pl, i_tr) = 0;
            if ris.tot.s_logico(i_pl, i_tr)
                % i seguenti parametri sono necessari per il calcolo
                % dell'armatura continua su tutta la lunghezza di L.Y
                l_.tr.anc.(supinf{i_si}) = 60*ris.tr.(['arm_' supinf{i_si}])(i_pl, i_tr).fi1 *1e-3;   % [m]lunghezza di ancoragggio per il diametro fi1 nell'ipotesi scarsa aderenza
                l_.tr.tra.(supinf{i_si})= 0.9*((2-i_si)*geom.tr.h + (2*(i_si==2) - 1) * arm.tr.(supinf{i_si}).d) * 1e-3; % [m] lunghezza di traslazione del momento
                l_.tr.beff.(supinf{i_si}) = L.b_max - l_.tr.anc.(supinf{i_si}) - l_.tr.tra.(supinf{i_si});   % [m] lunghezza effettiva delle barre
                l_.tr.Yeff.(supinf{i_si}) = L.Y/l_.tr.beff.(supinf{i_si}) * L.b_max;   % [m] lunghezza totale delle barre > di L.X perché considera anche le sovrapposizioni
                % calcola il peso delle barre fi1 i.e. l'armatura di base
                % continua. Il campo "peso_fi1" è esrpresso in [kg] pertanto è
                % necessario dividerlo per la lunghezza della campata per
                % ottenere l'incidenza al metro lineare che moltiplica la
                % lunghezza efficace l_.pl.Xeff.(...).
                ris.tr.peso.base(i_pl, i_tr) = ris.tr.peso.base(i_pl, i_tr) + l_.tr.Yeff.(supinf{i_si}) * ris.tr.(['arm_' supinf{i_si}])(i_pl, i_tr).peso_fi1 / ris.tr.(['arm_' supinf{i_si}])(i_pl, i_tr).L_fi1 * num.tr(i_pl) ;   % il peso è moltiplicato per L.Y/geom.pl.b per ottenere il totale esteso a tutta la superficie della platea
                % calcola il peso delle barre addizionali, non continue,
                % pertanto il peso è pari al peso_fi2
                ris.tr.peso.add(i_pl, i_tr) = ris.tr.peso.add(i_pl, i_tr) + ris.tr.(['arm_' supinf{i_si}])(i_pl, i_tr).peso_fi2 * (num.pali_tr(i_tr)-1);
                % Peso totale delle travi in funzione di L.tr
                ris.tr.peso.tot(i_pl, i_tr) = ris.tr.peso.base(i_pl, i_tr) + ris.tr.peso.add(i_pl, i_tr);
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
[C.tr, hw.tr] = contourf(L.pl, L.tr, ris.costo.tr');
clabel(C.tr, hw.tr)
title('Costo Travi')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Costo Travi')

subplot(3,2,5)
[C.pali, hw.pali] = contourf(L.pl, L.tr, ris.costo.pali');
clabel(C.pali, hw.pali);
title('Costo Pali')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Costo pali')

subplot(3,2,[2, 4, 6])
[C.tot, hw.tot] = contourf(L.pl, L.tr, ris.costo.tot'./costo_min, 1:.05:2);
clabel(C.tot, hw.tot);
title('Costo Totale [%]')
xlabel('L_{platea}')
ylabel('L_{trave}')
zlabel('Somma Costi')
%%
clc
fprintf('\nLa soluzione ottimale si ha in corripondenza degli indici:\n')
fprintf('\t\ti_min\t\tL[m]\n')
label = {'pl', 'tr'};
for i = 1:length(label)
    fprintf('%s\t% 8.0f\t% 8.2f \n', label{i}, i_min.(label{i}), L.(label{i})(i_min.(label{i})))
end
fprintf('\nChe corrispondono a %d travi e %d pali per trave, per un totale di %d pali.\n', num.tr(i_min.pl), num.pali_tr(i_min.tr), num.pali(i_min.pl, i_min.tr));
fprintf('\nLe deformazioni minime per questa combinazione sono pari a:\n');
fprintf('% -4s\t% 12s\t% 12s\t% s<%.2e\n', '', 's_min[mm]', 's_rel[%]', 's_rel', f_lim );
fprintf('% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'pl', ris.pl.s_min(i_min.pl), ris.pl.s_rel(i_min.pl), ris.pl.s_logico(i_min.pl));
fprintf('% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tr', ris.tr.s_min(i_min.pl, i_min.tr), ris.tr.s_rel(i_min.pl, i_min.tr), ris.tr.s_logico(i_min.pl, i_min.tr));
fprintf('% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot1', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel(i_min.pl, i_min.tr), ris.tot.s_logico(i_min.pl, i_min.tr));
fprintf('% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot2', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr)<= f_lim);