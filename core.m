%clearvars -except n n_sup n_inf risultati file2record
%close all
%clc
%% load
%load(file2record)
%% Aggiornamento template armatura
arm.pl.template = [arm.pl.cf, 0, 0; geom.pl.h - arm.pl.cf, 0, 0];  % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante

arm.tr.template = [arm.tr.cf, 0, 0; 280, 2, 14; 520, 2, 14; geom.tr.h - arm.tr.cf,  0, 0]; % [d, nb, fi] template della configurazione di armatura. La prima riga e l'ultima sono riscritte dalla funzione dimensionante


%% Variazione luce libera di inflessione della platea/trave
num.tr = num.campate + 1; % numeri di travi nella fondazione 10:29
L.pl = L.X ./ (num.tr - 1);      % [m] variazione della luce libera di inflessione della platea
L.tr = L.Y ./ (num.pali_tr - 1); % [m] variazione della luce libera di inflessione della trave

%% Dati dei pali
num.pali = num.tr' * num.pali_tr;   % numero totale di pali nelle varie configurazioni size: lenght(num.tr) x length(numpali_tr)
Fed.SLE_Q = L.pl' * L.tr * sollecitazioni.pl.q;
Fed.SLE_G = L.pl' * L.tr * (sollecitazioni.pl.q + 16.3) + ones(size(L.pl')) * (L.tr * (geom.tr.b * (geom.tr.h - geom.pl.h))*1e-6 * 25);
Fed.SLU = L.pl' * L.tr * sollecitazioni.pl.qSLU + ones(size(L.pl')) * (L.tr * 1.3 * (geom.tr.b * (geom.tr.h - geom.pl.h))*1e-6 * 25);

cedimento_globale.SLE_Q = Fed.SLE_Q ./ (Kw.ground * L.pl' * L.tr + Kw.palo);
cedimento_globale.SLE_G = Fed.SLE_G ./ (Kw.ground * L.pl' * L.tr + Kw.palo);
cedimento_globale.SLU = Fed.SLU ./ (Kw.ground * L.pl' * L.tr + Kw.palo);

RZ.palo.SLE_Q = Kw.palo * cedimento_globale.SLE_Q;
RZ.ground.SLE_Q = Kw.ground * L.pl' * L.tr .* cedimento_globale.SLE_Q;
RZ.ground.SLE_Qm2 = Kw.ground * cedimento_globale.SLE_Q;

RZ.palo.SLE_G = Kw.palo * cedimento_globale.SLE_G;
RZ.ground.SLE_G = Kw.ground * L.pl' * L.tr .* cedimento_globale.SLE_G;
RZ.ground.SLE_Gm2 = Kw.ground * cedimento_globale.SLE_G;

RZ.palo.SLU = Kw.palo * cedimento_globale.SLU;
RZ.ground.SLU = Kw.ground * L.pl' * L.tr .* cedimento_globale.SLU;
RZ.ground.SLUm2 = Kw.ground * cedimento_globale.SLU;

RZ.palo.cond = RZ.palo.SLE_G <= RZ.palo.max;

%% Inizio ciclo principale

% preallocamento variabili di carico
clear M
M.pl.inf = 0;
M.pl.sup = 0;
M.pl.mx = 0;
M.pl = repmat(M.pl, length(L.pl), 1);
M.tr.inf = 0;
M.tr.sup = 0;
M.tr.mx = 0;
M.tr = repmat(M.tr, length(L.pl), length(L.tr));
clear V
V.tr.mx = 0;
V.tr = repmat(V.tr, length(L.pl), length(L.tr));

%% SOLLECITAZIONI
clear soll
soll.pl.q = sollecitazioni.pl.q;
soll.pl.qSLU = sollecitazioni.pl.qSLU + 1.3 * geom.pl.peso_sezione;  % [kN/m2] carico agente in condizioni SLU
soll.pl.N = sollecitazioni.pl.N;
soll.pl.M = sollecitazioni.pl.M;   % funzione del momento sollecitante per trave 1 campata
soll.pl.V = sollecitazioni.pl.V; % taglio lungo l'asse della trave
% dati sollecitazione della trave
soll.tr.q = sollecitazioni.tr.q;  % inizializzo il valore di carico da stimare per ogni iterazione in funzione della luce di inflessione della platea
soll.tr.qSLU = sollecitazioni.tr.qSLU + 1.3 * sezione(geom.tr.b, geom.tr.h - geom.pl.h).peso_sezione;
soll.tr.N = sollecitazioni.tr.N;
soll.tr.M = sollecitazioni.tr.M;   % funzione del momento sollecitante per trave 1 campata
soll.tr.V = sollecitazioni.tr.V;
soll.tr = repmat(soll.tr, length(L.pl), 1);

% preallocamento risultati della freccia
clear ris
ris.pl.s_min = ones(size(length(L.pl)))*nan;
ris.pl.s_rel = ones(size(length(L.pl)))*nan;
ris.pl.s_logico = ones(size(length(L.pl)))*nan;
ris.pl.volume = geom.pl.area_sezione * 1e-6 * L.X * L.Y;


ris.tr.s_min = ones(length(L.pl), length(L.tr))*nan;
ris.tr.s_rel = ones(length(L.pl), length(L.tr))*nan;
ris.tr.s_logico = ones(length(L.pl), length(L.tr))*nan;
ris.tr.volume = num.tr * sezione(geom.tr.b, geom.tr.h - geom.pl.h).area_sezione * 1e-6 * L.Y;

ris.tot.s_min = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_rel = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_rel2 = ones(length(L.pl), length(L.tr))*nan;
ris.tot.s_logico = false(length(L.pl), length(L.tr));


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
clear armaturaValida
armaturaValida.pl = sign(sum(RZ.palo.cond, 2));
% armaturaValida.pl = repmat(armaturaValida.pl, length(L.pl), 1);
armaturaValida.tr = RZ.palo.cond;

% preallocamento risultati di costo
ris.costo.pl.steel = ones(length(L.pl),1) * nan;
ris.costo.pl.cls = ris.pl.volume * cu.cls;
ris.costo.pl.pompaggio = ris.pl.volume * cu.pompaggio;
ris.costo.pl.tot = ones(length(L.pl),1) * nan;

ris.costo.tr.steel = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.tr.cls = ris.tr.volume * cu.cls;
ris.costo.tr.pompaggio = ris.tr.volume * cu.pompaggio;
ris.costo.tr.casseri = 2 * num.tr * (sezione(geom.tr.b, geom.tr.h - geom.pl.h).h + 100) * 1e-3 * L.Y;
ris.costo.tr.tot = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.pali = ones(length(L.pl), length(L.tr)) * nan;
ris.costo.tot = ones(length(L.pl), length(L.tr)) * nan;


%% Calcolo della freccia
% Inizializzazione waitbar
clear ite
ite.current = 0;
ite.pl = 0;
ite.plMax = length(L.pl);
ite.trMax = length(L.tr);
ite.numMax = length(L.pl)*length(L.tr);
ite.progress = @(i) i/ite.numMax;
wb = waitbar(ite.progress(ite.current), sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));

% %% caricamento argomenti
% if exist('data_add.mat','file')
%     load('data_add.mat')
% end


%% corpo principale
for i_pl = 1:length(L.pl)
    
    %aggiornamento waitbar
    ite.pl = ite.pl + 1;
    ite.tr = 0;
    
    if armaturaValida.pl(i_pl)
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
        
        M.pl(i_pl).inf = k_Minf.pl * soll.pl.M(soll.pl.qSLU, L.pl(i_pl), L.pl(i_pl)/2);    % momento in mezzeria, in funzione della luce
        M.pl(i_pl).sup = k_Msup.pl * soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0);    % momento all'appoggio, in funzione della luce
        M.pl(i_pl).mx = soll.pl.M(soll.pl.qSLU, L.pl(i_pl), 0:dx:L.pl(i_pl));
        M.pl(i_pl).mx(M.pl(i_pl).mx > 0) = k_Minf.pl * M.pl(i_pl).mx(M.pl(i_pl).mx > 0);
        M.pl(i_pl).mx(M.pl(i_pl).mx < 0) = k_Msup.pl * M.pl(i_pl).mx(M.pl(i_pl).mx < 0);
        
        % dimensionamento dell'armatura superiore/inferiore
        [armaturaValida.pl(i_pl), ris_sup, ris_inf] = optimizeBeam(arm.pl, geom.pl.section, fck, M.pl(i_pl), soll.pl.N, L.pl(i_pl), opt1, opt2);
        if armaturaValida.pl(i_pl)
            ris.pl.arm_sup(i_pl) = ris_sup;
            ris.pl.arm_inf(i_pl) = ris_inf;
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
        % oltre al peso derivante dalla platea, si somma anche il peso proprio
        % della trave in c.a. (25.0 kN/m3), detratta la parte di sezione
        % sovrapposta alla platea (e.g. 0.5*(0.8-0.3)*25 )
        if armaturaValida.pl(i_pl)
            soll.tr(i_pl).q = soll.tr(i_pl).q + soll.pl.q * L.pl(i_pl);
            soll.tr(i_pl).qSLU = soll.tr(i_pl).qSLU + soll.pl.qSLU * L.pl(i_pl);   % carico agente sulle travi (carico x area di influenza)
        else
            soll.tr(i_pl).q = nan;
            soll.tr(i_pl).qSLU = nan;
        end
        
        % loop principale per i_tr
        
        for i_tr = 1:length(L.tr)
            
            % debug
            %             if i_pl == 5 && i_tr == 4
            %                 disp('yo')
            %             end
            
            %aggiornamento waitbar
            ite.tr = ite.tr + 1;
            ite.current = (ite.pl-1)*ite.trMax + ite.tr;
            waitbar(ite.progress(ite.current), wb, sprintf('Iterazione %d di %d\n%0.2f%%', ite.current, ite.numMax, ite.progress(ite.current)*100));
            set(wb, 'Name', sprintf('Progress \n%0.2f%%', ite.progress(ite.current)*100));
            
            if armaturaValida.tr(i_pl, i_tr)
                
                % sollecitazione agente
                M.tr(i_pl, i_tr).inf = k_Minf.tr * soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), L.tr(i_tr)/2);    % momento in mezzeria, in funzione della luce
                M.tr(i_pl, i_tr).sup = k_Msup.tr * soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0);    % momento all'appoggio, in funzione della luce
                M.tr(i_pl, i_tr).mx = soll.tr(i_pl).M(soll.tr(i_pl).qSLU, L.tr(i_tr), 0:dx:L.tr(i_tr));
                M.tr(i_pl, i_tr).mx(M.tr(i_pl, i_tr).mx > 0) = k_Minf.tr * M.tr(i_pl, i_tr).mx(M.tr(i_pl, i_tr).mx > 0);
                M.tr(i_pl, i_tr).mx(M.tr(i_pl, i_tr).mx < 0) = k_Msup.tr * M.tr(i_pl, i_tr).mx(M.tr(i_pl, i_tr).mx < 0);        
                V.tr(i_pl, i_tr).vx = soll.tr(i_pl).V(soll.tr(i_pl).qSLU, L.tr(i_tr), 0:dx:L.tr(i_tr));
                
                % dimensionamento dell'armatura superiore/inferiore
                [armaturaValida.tr(i_pl, i_tr), ris_sup, ris_inf] = optimizeBeam(arm.tr, geom.tr.section, fck, M.tr(i_pl, i_tr), soll.tr(i_pl).N, L.tr(i_tr), opt1, opt2);
                if armaturaValida.tr(i_pl, i_tr)
                    ris.tr.arm_sup(i_pl, i_tr) = ris_sup;
                    ris.tr.arm_inf(i_pl, i_tr) = ris_inf;
                end
                
                % calcolo dell'armatura a taglio
                if armaturaValida.tr(i_pl, i_tr)
                    [armSt_, pesoSt_]  = calcPesoStaffe(geom.tr.section, L.tr(i_tr)*1e3, arm.tr.template(:,1), arm.tr.staffe.nb_sw, arm.tr.staffe.fi_sw, min(ris_inf.fi1, ris_sup.fi1), arm.tr.staffe.s_lim, fck, soll.tr(i_pl).N, V.tr(i_pl, i_tr).vx, 'copriferro', arm.tr.staffe.cf);
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
end
close(wb)
clear('h')

