function [ ris ] = dimenTaglio( sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, Ved, varargin)
%DIMENTAGLIO Dimensionamento della trave a taglio
%   calcolo delle differenti configurazioni di armatura a taglio per la
%   trave di lunghezza L, soggetta al taglio Ved. In funzione delle
%   combinazioni di armatura a taglio (nb, fi, s) si ricavano i tagli
%   resistenti t.c. Vrd >= Ved.
%   Argomenti di input:
%       sezione:    matrice Nx4 [xm, ym, dx, dy] che discretizza la sezione
%       d:  altezza utile della sezione in c.a.
%       nb_sw: numero di bracci (vettore delle diverse possibilità)
%       fi_sw: vettore dei diametri dell'armatura trasversale
%       fi_sl: diametro delle barre longitudinali
%       s_lim: variabile del passo. Se l'opzione passo == 'variabile'
%       (default) allora s_lim(1) è il passo minimo possibile, e s_lim(2) è
%       la variazione del passo ad ogni iterazione (ds)
%       ammissibile, ed il secondo elemento è la risoluzione del passo
%       fck: resistenza caratteristica del cls
%       Ved: taglio sollecitante [kN]
%
%   Argomenti opzionali:
%       alpha: [90] inclinazione dell'armatura a taglio
%       precisione: [12] precisione del calcolo
%       passo: ['variabile'] opzione che determina l'utilizzo del parametro
%       s_lim
%       zona_critica: [true] condizione per il calcolo del passo massimo

%% Valori di default ed estrazione argomenti opzionali
%lock = false;
alpha = 90; % inclinazione dell'armatura a taglio
classe = 'B';
zona_critica = true;
passo = 'variabile';

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'alpha'
            alpha = varargin{2};
        case 'classe'
            classe = varargin{2};
        case 'zona_critica'
            zona_critica = varargin{2};
        case 'passo'
            passo = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% estrazione dei vettori della matrice "sezione"

xm = sezione(:,1);  % il primo campo alla coordinata in xm del baricentro dei rettangoli
ym = sezione(:,2);  % il secondo campo è riservato alla coordinata in ym del baricentro dei rattngoli
db = sezione(:,3);  % il terzo campo è riservato alla larghezza in x dei rettangoli
dh = sezione(:,4);  % il quarto campo è riservato alla larghezza in y dei rettangoli

[uym, ih, ~] = unique(ym);    % elimina le componenti duplicate
H = sum(dh(ih));    % massima altezza della sezione (somma solo una volta dh alla quota ym)
B = zeros(length(uym),2);
for i = 1:length(uym)
    B(i,:) = [sum( abs( db( ym == uym(i) ) ) ), uym(i)]; % larghezza della sezione in funzione dell'altezza
end

bw = min(B(:,1));   % la larghezza della sezione resistente a taglio è il minimo delle larghezze della sezione

%% funzioni locali
As.fun = @(n,fi) n*pi*fi^2/4; % area di armatura per il diametro e numero di barre
func = ['verificaTaglio(sezione, dmax, Asw, ~, s, alpha, fck, Ned, ''classe'', ''' classe ''');'];
dmax = max(d);  % altezza utile della sezione


%% corpo principale
i_r = 0;    % iteratore per la struttura dei risultati
ris = struct('nb', nan, 'fi', nan, 'Asw', nan, 's', nan, 's_max', nan, 'Asw_s', nan,...
    'cot_theta', nan, 'Vrsd', nan, 'Vrcd', nan, 'Vrd', nan, 'ratio', nan);
ris = repmat(ris, length(nb_sw)*length(fi_sw),1);

for i_nb = 1:length(nb_sw)
    soluzione_finale = false;
    for i_fi = 1:length(fi_sw)
        Asw = As.fun(nb_sw(i_nb), fi_sw(i_fi));
        if zona_critica
            s_max = min([dmax/4, 225, 8*fi_sl, 24*fi_sw(i_fi)]); % passo massimo in zona critica
        else
            s1 = Asw/(1.5*bw)*1e3;  % passo in funzione dell'armatura a taglio
            s_max = min([s1, 1e3/3, 0.8*dmax]); % passo massimo in zona non critica
        end
        
        % passo
        if strcmp(passo, 'variabile')
            s_min = s_lim(1);
            ds = s_lim(2);
            s = floor(s_max/ds)*ds;
            j = 0;  % iteratore di ciclo per passo variabile
        else
            indice_s_validi = s_lim <= s_max;
            s_ = sort(s_lim(indice_s_validi), 'descend');    % riduce il vettore s a solo i rultati validi
            s_min = s_(end);
        end
        Vrd = 0; % variabile di controllo del ciclo
        while Vrd < abs(Ved)
            if strcmp(passo, 'variabile')
                s = s - ds*j;   % il passo s viene diminuito ad ogni ciclo del valore di ds
            else
                s = s_(1);  % viene assegnato il primo valore della variabile temporanea ad s
            end
            % condizione di fine ciclo
            if s <= s_min || s <= 0
                break
            end
            [Vrd, Vrsd, Vrcd, cot_theta] = verificaTaglio(sezione, dmax, Asw, 0, s, alpha, fck, Ned, 'classe', classe);
            Vrd = Vrd*1e-3;
            Vrsd = Vrsd*1e-3;
            Vrcd = Vrcd*1e-3;
            
            % salvataggio dei dati per la soluzione corretta
            if Vrd >= abs(Ved)
                i_r = i_r+1;    % aggiornamento iteratore
                ris(i_r).nb = nb_sw(i_nb);
                ris(i_r).fi = fi_sw(i_fi);
                ris(i_r).Asw = Asw;
                ris(i_r).s = s;
                ris(i_r).s_max = s_max; % passo massimo per Asw corrente
                ris(i_r).Asw_s = Asw/s;
                ris(i_r).cot_theta = cot_theta;
                ris(i_r).Vrsd = Vrsd;
                ris(i_r).Vrcd = Vrcd;
                ris(i_r).Vrd = Vrd;
                ris(i_r).ratio = abs(Ved)/Vrd;
            end
            
            % aggiornamento delle varibili di ciclo
            if strcmp(passo, 'variabile')
                j = j+1;
                
                % condizione di fine ciclo
                if Vrd >= abs(Ved) && s == floor(s_max/ds)*ds
                    % se s è un valore accettabile (Vrd >= Ved) e s è apri
                    % al più grande valore ammissibile, allora la soluzione
                    % trovata è la migliore per la combinazione (fi,s). Non
                    % ha senso cercare altre soluzione perché si avrebbe un
                    % diametro maggiore a passo costante, dunque soluzione
                    % meno efficiente.
                    soluzione_finale = true;
                end
            else
                s_(1) = [];
                if isempty(s_)
                    break
                end
                
                % condizione di fine ciclo
                if Vrd >= abs(Ved) && s == s(1) 
                    % se s è un valore accettabile (Vrd >= Ved) e s è apri
                    % al più grande valore ammissibile, allora la soluzione
                    % trovata è la migliore per la combinazione (fi,s). Non
                    % ha senso cercare altre soluzione perché si avrebbe un
                    % diametro maggiore a passo costante, dunque soluzione
                    % meno efficiente.
                    soluzione_finale = true;
                end
            end
        end
        if soluzione_finale
            break
        end
    end
end
%% eliminazione delle righe vuote della struct
while isnan(ris(end).nb)
        ris(end) = [];
end


end

