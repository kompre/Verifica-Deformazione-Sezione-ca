%% ANALISI DEI RISULTATI
% Combina i risultati delle frecce derivanti da core e calcola i costi
% delle diverse soluzioni.

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
            ris.tot.s_logico(i_pl, i_tr) = ris.tot.s_rel(i_pl, i_tr) <= f_lim && ris.pl.s_rel(i_pl) <= f_lim && ris.tr.s_rel(i_pl, i_tr) <= f_lim;
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
        ris.costo.pl.steel(i_pl) = ris.pl.peso.tot(i_pl) * cu.steel;
        ris.costo.pl.tot(i_pl) = ris.costo.pl.steel(i_pl) + ris.costo.pl.cls + ris.costo.pl.pompaggio;
        for i_tr = 1:length(L.tr)
            if ris.tot.s_logico(i_pl, i_tr)
                ris.costo.tr.steel(i_pl, i_tr) = ris.tr.peso.tot(i_pl, i_tr) * cu.steel;
                ris.costo.tr.tot(i_pl, i_tr) = ris.costo.tr.steel(i_pl, i_tr) + ris.costo.tr.cls(i_pl) + ris.costo.tr.pompaggio(i_pl) + ris.costo.tr.casseri(i_pl);
                ris.costo.pali(i_pl, i_tr) = num.pali(i_pl, i_tr) * cu.palo;
                ris.costo.tot(i_pl, i_tr) = ris.costo.pl.tot(i_pl) + ris.costo.tr.tot(i_pl, i_tr) + ris.costo.pali(i_pl, i_tr);
            end
        end
    end
end

[costo_min, I] = min(ris.costo.tot(:));
[i_min.pl, i_min.tr] = ind2sub(size(ris.costo.tot), I);
ris.costo.min = costo_min;
ris.costo.ind_min = [i_min.pl, i_min.tr];