function [ armaturaValida, varargout ] = optimizeBeam( arm, sezione, fck, Med, Ned, L, opt1, opt2 )
%OTTIMIZZAZIONECAMPATA dimensiona l'armatura a flessione della trave di
%lunghezza L sollecitata da Med, Ned.
%  Dimensiona e ottimizza la sezione di c.a. sia a momento positivo (inf)
%  che a momento negativo (sup). L'ottimo è ricavato tramite un processo di
%  iterazione.
%  Parametri:
%   arm: matrice Nx[d, nb1, fi1] che specifica le diverse
%       configurazioni di armatura in funzione di d; La prima riga sarà
%       riscritta dal dimensionamento dell'armatura superiore, mentre l'ultima
%       riga da quella dell'armatura inferiore.
%   sezione: matrice Nx[xm, ym, dx, dy] che discretizza la sezione di c.a.
%       in esame.
%   fiv: struttura di vettori che contiene i diametri ammissibili. (.inf/.sup)
%   nb1: struttura di vettori che contiene la quantità di barre per riga (.inf/.sup)
%   ammissibili per fi1. nb2: vettore che contiene la quantità di barre per
%   riga ammissibili per fi2. fck: resistenza caratteristica del cls. Med:
%   struttura che contiene i campi 'inf', 'sup', e 'mx'. Ned: sforzo
%   normale. opt1: opzioni da fornire alla funzione dimenSezione(...,
%   varargin). opt2: opzioni da fornire alla funzione calcPesoCamp(...,
%   varargin)


if ~sum(strcmp(opt1,'lock'))
    opt1 = [opt1, {'lock', 'false'}];
end
indLock = find((strcmp(opt1,'lock'))); % determina l'indice dove si trova lock

% assumo che l'armatura sia valida
armaturaValida = true;

varargout = cell(1,2);


% dimensionamento dell'armatura superiore/inferiore
if abs(Med.inf) <= abs(Med.sup)
    lembo = {'inf', 'sup'};
    l_sos = [length(arm.template(:,1)), 1];
else
    lembo = {'sup', 'inf'};
    l_sos = [1, length(arm.template(:,1))];
end

j = 0; % variabile di ciclo while
continua_ciclo = true;
while continua_ciclo % continua_ciclo.sup == true || continua_ciclo.inf == true
    i_lem = 1 + rem(j,2);
    
    % aggiornamento dell'armatura tipo
    if j == 0
        arm.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [0, 0];    % per la prima iterazione, l'armatura è semplice
    else
        arm.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [ris_.nb1, ris_.fi1];    % aggiorno con i valori dell'iterazione precedente
    end
    
    opt1{indLock+1} = arm.(lembo{i_lem}).lock; % aggiorna il valore di lock nelle opzioni secondi il lembo corrente
    % calcolo momento resistente e peso
    ris_ = dimenSezione(sezione, arm.template, arm.(lembo{i_lem}).fiv, arm.(lembo{i_lem}).nb1, arm.(lembo{i_lem}).nb2, fck, Med.(lembo{i_lem}), Ned, opt1{:});
    ris_ = calcPesoCamp(ris_, sign(Med.(lembo{i_lem})), L, Med.mx, opt2{:});
    
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
    lembo_opposto = lembo{~strcmp(lembo{i_lem},lembo)};
    if   j >= 2 && ~isempty(ris_) && ris_.As_tot == ris.(['arm_' lembo{i_lem}]).As_tot
        % condizione di fine ciclo: l'armatura è rimasta invariata dal
        % ciclo precedente
        continua_ciclo = false;
    elseif j >= 2 && ~isempty(ris_) && ris_.As_tot < ris.(['arm_' lembo{i_lem}]).As_tot
        % l'armatura è diminuita rispetto al ciclo precedente, quindi
        % si continua l'iterazione. Inoltre salvo l'iterazione
        % precedente in memoria
        continua_ciclo = true;
    elseif j >= 2 && (isempty(ris_) || ris_.As_tot > ris.(['arm_' lembo{i_lem}]).As_tot && exist('ris2', 'var') && isfield(ris2, ['arm_' lembo{i_lem}]))
        % l'armatura è aumentata rispetto al ciclo precedente.
        % L'armatura calcolata per il ciclo precedente è la migliore,
        % pertanto salvo quella soluzione ed interrompo il ciclo.
        continua_ciclo = false;
        
        % ripristino l'armatura calcolata due iterazioni prima per
        % il lembo opposto, e calcolo l'armatura per il lembo
        % attuale imponendo la configurazione di armatura
        opt3 = {'tipo', '''elastica''', 'lock', 'false', 'precisione', '6', 'ottimizzazione', 'false'};
        if strcmp(lembo{1}, lembo{i_lem})
            fiv_ = [ris2.(['arm_' lembo{i_lem}]).fi1, ris2.(['arm_' lembo{i_lem}]).fi2];
            nb1_ = ris2.(['arm_' lembo{i_lem}]).nb1;
            nb2_ = ris2.(['arm_' lembo{i_lem}]).nb2;
            arm.template(l_sos(~strcmp(lembo{i_lem},lembo)),2:3) = [ris2.(['arm_' lembo_opposto]).nb1, ris2.(['arm_' lembo_opposto]).fi1];    % aggiorno con i valori dell'iterazione precedente
            ris.(['arm_' lembo_opposto]) = ris2.(['arm_' lembo_opposto]);
            ris_ = dimenSezione(sezione, arm.template, fiv_, nb1_, nb2_, fck, Med.(lembo{i_lem}), Ned, opt3{:});
            ris_ = calcPesoCamp(ris_, sign(Med.(lembo{i_lem})), L, Med.mx, opt2{:});
        else
            fiv_ = [ris2.(['arm_' lembo_opposto]).fi1, ris2.(['arm_' lembo_opposto]).fi2];
            nb1_ = ris2.(['arm_' lembo_opposto]).nb1;
            nb2_ = ris2.(['arm_' lembo_opposto]).nb2;
            arm.template(l_sos(~strcmp(lembo_opposto,lembo)),2:3) = [ris2.(['arm_' lembo{i_lem}]).nb1, ris2.(['arm_' lembo{i_lem}]).fi1];    % aggiorno con i valori dell'iterazione precedente
            ris_ = dimenSezione(sezione, arm.template, fiv_, nb1_, nb2_, fck, Med.(lembo_opposto), Ned, opt3{:});
            ris_ = calcPesoCamp(ris_, sign(Med.(lembo_opposto)), L, Med.mx, opt2{:});
            ris.(['arm_' lembo_opposto]) = ris_;
            ris_ = ris.(['arm_' lembo{i_lem}]);
        end
    end
    
    if ~isempty(ris_)
        if j >= 2
            ris2.(['arm_' lembo{i_lem}]) = ris.(['arm_' lembo{i_lem}]);
        end
        ris.(['arm_' lembo{i_lem}]) = ris_; % salvo la soluzione temporanea nel struttura globale
    else
        armaturaValida = false; % flag che segnala se la combinazione di armatura per la luce in oggetto è valida oppure no
        break
    end
    j = j+1;
end
if armaturaValida
    varargout{1} = ris.arm_sup;
    varargout{2} = ris.arm_inf;
end
end

