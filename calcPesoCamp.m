function [ris] = calcPesoCamp( ris, segno, L, Med, varargin)
%CALCOLOPESO calcola il peso delle barre di acciaio sulla trave  di
% lunghezza L per la configurazione di armatura contenuta in ris.
% Il calcolo è effettuato nelle ipotesi che sia previsti due tipi di
% armatura: quella base (1) che è applicata lunga tutta la lunghezza della
% trave L, e quella aggiuntiva (2), applicata solo dove è necessaria, i.e.
% Med > Mrd_fi1
% Per semplificare il processo Med è un vettore dei valori del momento
% sollecitante, lungo x
%   ris è una struttura che contiene i seguenti campi:
%       nb1: numero di barre di diametro 1;
%       nb2: numero di barre di diametro 2;
%       fi1: diametro delle barre di base;
%       fi2: diametro delle barre aggiuntive;
%       As1: area armatura di base;
%       As2: area armatura aggiuntiva;
%       Mrd_fi1: momento resistente base;
%       Mrd:    momento resistente armatura totale;
%       ratio: rapporto tra Med/Mrd;
%       armatura: struttura che riassume tutti i livelli dell'armatura
%       previsti in funzione dell'altezza utile d

%% Valori di default ed estrazione argomenti opzionali
peso_acciaio = 7850;
% parametri per lunghezza di ancoraggio
a1 = 1;     % barre non piegate
a2 = 60;    % scarsa aderenza

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'aderenza'
            if strcmp(varargin{2},'buona')
                a2 = 40;
            end
        case 'piegate'
            if varargin{2}
                a1 = 0.7;
            end
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% Calcolo della lunghezza delle barre per armatura aggiuntiva (nb2, fi2)

if isempty(ris)
    return
end


dx = L/(length(Med)-1);
x = 0:dx:L;
for i_r = 1:length(ris)
    % se segno < 0 allora si sta verificando una sezione negativa per cui
    % il momento sollecitante è maggiore di quello resistente quando
    % Med<Mrd (entrambi negativi).
    % se segno > 0 allora si sta verificando una sezione positiva quindi,
    % il momento sollecitante è maggiore di quello resistente quando
    % Med>mrd
    if segno < 0
        cond = Med < segno * ris(i_r).Mrd_fi1;
    else
        cond = Med > segno * ris(i_r).Mrd_fi1;
    end
    
    j = 0;  % variabile di ciclo while
    while ~isempty(cond);
        i_fineConcio = find(cond == ~cond(1), 1, 'first')-1;    % trova il primo valore diverso dal primo elemento (e sottrae 1, i.e. torva l'ultimo elemento in cui tutti i valori sono 1 o 0)
        if isempty(i_fineConcio)
            i_fineConcio = length(cond);
        end
        j = j+1;
        if j == 1
            conci(j,:) = [1, i_fineConcio, cond(1)];   % primo concio
        else
            conci(j,:) = [conci(j-1,2)+1, conci(j-1,2)+i_fineConcio, cond(1)];   % memorizza gli estremi dei conci
        end
        cond(1:i_fineConcio) = [];  % elimina il concio trovato in questo ciclo
    end
    conci = conci( conci(:,3)==1,:);    % elimino i conci dove Mrd>Med
    
    % Lunghezza di ancoraggio
    % a1: fattore di riduzione per barre piegate
    % a2: moltiplicatore del diametro in funzione delle condizioni di
    % aderenza
    la = a1 * a2 * ris(i_r).fi2 *1e-3;
    
    % lunghezza di traslazione dei momenti
    lt = 0.9*max(ris(i_r).armatura.d)/2 *1e-3;
    
    % estrazione delle x
    [row,~] = size(conci);
    lx = zeros(row, 1);
    L_tot = 0;
    ris(i_r).Lp_fi2 = struct;
    for j = 1:row
        % determina se il concio è posto alle estremità della trave, oppure
        % se è interno. Nel caso sia alle estremità, la lunghezza di
        % ancoraggio si conta 1 sola volta (n=1) se invece è posto
        % internamente, allora la lunghezza di ancoraggio si conta due
        % volta  (n==2)
        if or(conci(j,1) == 1, conci(j,2) == length(x))
            n = 1;
        else
            n = 2;
        end
        lx(j) = x(conci(j,2)) - x(conci(j,1)) + n*(la + lt); % lunghezza del tratto da coprire con armatura addizionale tale per cui Med = Mrd_base
        ris(i_r).Lp_fi2.(['Lp' num2str(j)]) = lx(j);    % salva le lunghezze parziali nella struttura dei risultati
        L_tot = L_tot + lx(j);  % somma di tutte le lunghezze parziali
    end
    ris(i_r).L_fi1 = L;
    ris(i_r).L_fi2 = L_tot;
    ris(i_r).peso_fi1 = ris(i_r).As1*1e-6 * ris(i_r).L_fi1 * peso_acciaio;
    ris(i_r).peso_fi2 = ris(i_r).As2*1e-6 * ris(i_r).L_fi2 * peso_acciaio;
    ris(i_r).peso_tot = ris(i_r).peso_fi1+ris(i_r).peso_fi2;
end

%% calcolo del peso per intera sezione

for i_r = 1:length(ris)
    % assegnazione dei valori nuovi
    ris(i_r).armatura.la_fi1 = 1 * 60 * ris(i_r).armatura.fi1*1e-3;  % lunghezza di ancoraggio in scarsa aderenza. Siccome è per le barre continue, a1 = 1, i.e. non piegate
    ris(i_r).armatura.lt_fi1(:) = sign(ris(i_r).armatura.nb1) .* lt;  % traslazione dei momenti
    ris(i_r).armatura.L_fi1 = sign(ris(i_r).armatura.nb1) .* ris(i_r).L_fi1;
    ris(i_r).armatura.L_fi2 = sign(ris(i_r).armatura.nb2) .* ris(i_r).L_fi2;
    ris(i_r).armatura.peso_fi1 = ris(i_r).armatura.As1*1e-6 .* ris(i_r).armatura.L_fi1 * peso_acciaio;
    ris(i_r).armatura.peso_fi2 = ris(i_r).armatura.As2*1e-6 .* ris(i_r).L_fi2 * peso_acciaio;
    ris(i_r).armatura.peso_tot = ris(i_r).armatura.peso_fi1 + ris(i_r).armatura.peso_fi2;
    ris(i_r).armatura.Peso = sum( ris(i_r).armatura.peso_tot );
end


