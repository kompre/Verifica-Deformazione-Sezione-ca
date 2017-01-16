function [ arm, peso_totale] = calcPesoStaffe( sezione, L, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, vx, varargin)
%CALCPESOSTAFFE calcolo accurato del peso delle staffe per trave di luce L
%   Calcola le staffe necessarie per una trave di luce L, soggetta a taglio
%   sollecitante vx.
%   Parametri di input:
%       vx: vettore dei valori di taglio in funzione di x

%% Valori di default ed estrazione argomenti opzionali
alpha = 90; % inclinazione dell'armatura a taglio
classe = 'B';
passo = 'variabile';
copriferro = 25;

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'alpha'
            alpha = varargin{2};
        case 'classe'
            classe = varargin{2};
        case 'passo'
            passo = varargin{2};
        case 'copriferro'
            copriferro = varargin{2};
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

%%

x = linspace(0, L, length(vx)); % vettore delle ascisse associato a vx
x_cr_sx = H; % larghezza della zona critica dall'appoggio sinistro
x_cr_dx = L-H; % larghezza della zona critica dall'appoggio destro
indSt.zcr.sx = x <= x_cr_sx;
indSt.zcr.dx = x >= x_cr_dx;
indSt.zcr.comb = and(indSt.zcr.sx, indSt.zcr.dx);
indTr.sx = not([ones(size(x(1:floor(length(x)/2))))*false, ones(size(x(floor(length(x)/2)+1:end)))*true]); % indice logico che individua la metà sinistra della trave
indTr.dx = not(indTr.sx); % indice logico che individua la metà destra della trave

% determinazione dei valori di taglio fondamentali
V.sx = abs(vx(1));   % taglio all'appoggio sinistro
V.dx = abs(vx(end)); % taglio all'appoggio destro

V.cr_sx = min(abs(vx(indSt.zcr.sx)));  % taglio al limite della zona critica a sx
V.cr_dx = min(abs(vx(indSt.zcr.dx)));  % taglio al limite della zona critica a sx

% dimensionamento armatura per la zona critica
lato = {'sx','dx'};
for i_lato = 1:length(lato)
    % dimensiono l'armatura a taglio e salvo la stuttura Nx1 in una
    % variabile temporanea
    try
        arm_ = dimenTaglio(sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, V.(lato{i_lato}), 'classe', classe, 'zona_critica', true, 'passo', passo);
    catch
        arm = [];   % se il dimensionamento fallisce, allora torna una struct vuota
        peso_totale = nan;
        return
    end
    
    % ricava armatura ottimale (criterio minore Asw/s)
    while length(arm_) > 1
        if arm_(1).Asw_s < arm_(2).Asw_s
            arm_(2) = [];   % elimino il risultato maggiore
        else
            arm_(1) = [];   % elimino il risultato maggiore
        end
    end
    zcr.(lato{i_lato}) = arm_;   % salvo il minimo
    
    % interrompo il ciclo se V.sx == V.dx
    if V.sx == V.dx     % se vera, questa condizione interrompe il ciclo alla prima iterazione. Se falsa, lo sarà sempre
        zcr.dx = zcr.sx;
        break
    end
end



for i_lato = 1:length(lato)
    % dimensiono l'armatura a taglio e salvo la stuttura Nx1 in una
    % variabile temporanea
    try
        arm_ = dimenTaglio(sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, V.(['cr_', lato{i_lato}]), 'classe', classe, 'zona_critica', false, 'passo', passo);
    catch
        return
    end
    
    % ricava armatura ottimale (criterio minore Asw/s)
    while length(arm_) > 1
        if arm_(1).Asw_s < arm_(2).Asw_s
            arm_(2) = [];   % elimino il risultato maggiore
        else
            arm_(1) = [];   % elimino il risultato maggiore
        end
    end
    znc.(lato{i_lato}) = arm_;   % salvo il minimo
    
    % ricavo il valore del taglio in funzione di Asw/s_max
    indStBase.(lato{i_lato}) = ones(size(x))*false; % inzializzo l'indice delle staffe di base
    
    % converto s_max in un valore reale, in funzione del tipo di passo.
    if strcmp(passo, 'variabile')
        smax_ = floor(arm_.s_max/s_lim(2))*s_lim(2);
    else
        smax_ = s_lim(end);
    end
    
    % Calcolo del taglio resistente per smax_
    if arm_.s ~= smax_
        [Vrd_smax, Vrsd_smax, Vrcd_smax, cot_theta_smax] = verificaTaglio(sezione, d, arm_.Asw, 0, smax_, alpha, fck, Ned, 'classe', classe);
         
        % salvo l'armatura di base, che ha stesso numero di braccia e diametro di quella minima
        armBase.(lato{i_lato}) = arm_;
        
        % aggiorno i campi necessari (attenzione: la funzione
        % 'verificaTaglio' ritorna risultati coerenti con le dimensioni dei
        % dati inseriti [MPa, mm] --> [N], pertanto è necessario scalare il
        % risultato per ottenere valori congruenti. La funzione dimenTaglio
        % realizza già questa conversione.
        armBase.(lato{i_lato}).s = smax_;
        armBase.(lato{i_lato}).Asw_s = armBase.(lato{i_lato}).Asw/armBase.(lato{i_lato}).s;
        armBase.(lato{i_lato}).Vrd = Vrd_smax*1e-3;
        armBase.(lato{i_lato}).Vrsd = Vrsd_smax*1e-3;
        armBase.(lato{i_lato}).Vrcd = Vrcd_smax*1e-3;
        armBase.(lato{i_lato}).cot_theta = cot_theta_smax;
    else
        armBase.(lato{i_lato}) = arm_;    % se s == s_max, allora il taglio minimo è pari al taglio resistente di progetto
    end
    
    % indice logico delle ascisse dove Ved <= Vrd_base
    indStBase.(lato{i_lato}) = abs(vx) <= armBase.(lato{i_lato}).Vrd;   % questo indice contiene anche valori già coperti dall'armatura della zona critica
    indStBase.(lato{i_lato}) = and(and(~indSt.zcr.(lato{i_lato}), indStBase.(lato{i_lato})), indTr.(lato{i_lato})); % elimina i valori della zona critica ed isola solo quelli per il lato destro o sinistro
end

% combina i due indici in uno unico perché come armatura di base si prende
% la maggiore di quella determinata per il lato sinistro e quella
% determinata per il lato destro.
indSt.base.comb = or(indStBase.sx, indStBase.dx);    % combina i due indici in un unico indice
if armBase.sx.Vrd < armBase.dx.Vrd
    armBase.comb = armBase.dx;
else
    armBase.comb = armBase.sx;
end

% calcolo dell'indice per la zona non critica

indSt.znc.sx = ~or(or(indSt.zcr.sx, indSt.base.comb), indTr.dx);
indSt.znc.dx = ~or(or(indSt.zcr.dx, indSt.base.comb), indTr.sx);

indSt.tot = [sum(indSt.zcr.sx), sum(indSt.znc.sx), sum(indSt.base.comb), sum(indSt.znc.dx), sum(indSt.zcr.dx)] > 0;
indSt.comb = [indSt.zcr.sx; indSt.znc.sx; indSt.base.comb; indSt.znc.dx; indSt.zcr.dx];

arm = [zcr.sx, znc.sx, armBase.comb, znc.dx, zcr.sx];
for i = 1:length(indSt.tot)
    if indSt.tot(i)
        x_rid = x(indSt.comb(i,:));
        l = x_rid(end) - x_rid(1);
        if i == 1 || i == length(indSt.tot)
            arm(i).zona = 'critica';
        else
            arm(i).zona = 'non critica';
        end
        arm(i).L_zona = l;
        arm(i).Q_st = round(arm(i).L_zona/arm(i).s);
        [~, peso_] = lunghezzaStaffa(H, bw, arm(i).nb, copriferro, 'fi', sprintf('%d', arm(i).fi));
        arm(i).Peso = arm(i).Q_st * peso_;
    end
end

%% eliminazione righe vuote
j = 1;
while j <= length(arm)
    if isempty(arm(j).zona)
        arm(j) = [];
    else
        j = j+1;
    end
end

%% calcolo del peso totale
peso_totale = 0;
for i = 1:length(arm)
    peso_totale = peso_totale + arm(i).Peso;
end

end




