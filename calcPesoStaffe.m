function [ arm ] = calcPesoStaffe( H, L, vx, sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, varargin)
%CALCPESOSTAFFE calcolo accurato del peso delle staffe per trave di luce L
%   Calcola le staffe necessarie per una trave di luce L, soggetta a taglio
%   sollecitante vx.
%   Parametri di input:
%       vx: vettore dei valori di taglio in funzione di x

%% Valori di default ed estrazione argomenti opzionali
peso_acciaio = 7850;
alpha = 90; % inclinazione dell'armatura a taglio
classe = 'B';
passo = 'variabile';

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'alpha'
            alpha = varargin{2};
        case 'classe'
            classe = varargin{2};
        case 'passo'
            passo = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

x = linspace(0, L, length(vx)); % vettore delle ascisse associato a vx
x_cr_sx = H; % larghezza della zona critica dall'appoggio sinistro
x_cr_dx = L-H; % larghezza della zona critica dall'appoggio destro
indSt.zcr.sx = x <= x_cr_sx;
indSt.zcr.dx = x >= x_cr_dx;
indSt.zcr.comb = and(indSt.zcr.sx, indSt.zcr.dx);
indTr.sx = not([ones(size(x(1:floor(length(x)/2))))*false, ones(size(x(ceil(length(x)/2):end)))*true]); % indice logico che individua la metà sinistra della trave
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
    arm_ = dimenTaglio(sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, V.(lato{i_lato}), 'classe', classe, 'zona_critica', true, 'passo', passo);
    
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
    arm_ = dimenTaglio(sezione, d, nb_sw, fi_sw, fi_sl, s_lim, fck, Ned, V.(['cr_', lato{i_lato}]), 'classe', classe, 'zona_critica', false, 'passo', passo);
    
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
    if arm_.s ~= arm_.s_max
        [Vrd_smax, Vrsd_smax, Vrcd_smax, cot_theta_smax] = verificaTaglio(sezione, d, arm_.Asw, 0, arm_.s_max, alpha, fck, Ned, 'classe', classe);
        indStBase.(lato{i_lato}) = abs(vx) <= Vrd_smax *1e-3;
        indStBase.(lato{i_lato}) = and(and(~indSt.zcr.(lato{i_lato}), indStBase.(lato{i_lato})), indTr.(lato{i_lato}));
        
        % salvo l'armatura di base, che ha stesso numero di braccia e diametro di quella minima
        armBase.(lato{i_lato}) = arm_;
        % aggiorno i campi necessari
        armBase.(lato{i_lato}).s = arm_.s_max;
        armBase.(lato{i_lato}).Vrd = Vrd_smax*1e-3;
        armBase.(lato{i_lato}).Vrsd = Vrsd_smax*1e-3;
        armBase.(lato{i_lato}).Vrcd = Vrcd_smax*1e-3;
        armBase.(lato{i_lato}).cot_theta = cot_theta_smax;
    else
        armBase.(lato{i_lato}) = arm_;    % se s == s_max, allora il taglio minimo è pari al taglio resistente di progetto
    end
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
for i = 1:length(indSt.tot)
    if indSt.tot(i)
        x_rid = x(indSt.comb(i,:));
        l = x_rid(end) - x_rid(1)
    end
end





