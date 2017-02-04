function [ ris ] = dimenSezione( sezione, armatura, fiv, Nb1, Nb2, fck, Med, Ned, varargin )
%DIMENBARRE dimensiona le barre in funzione del momento sollecitante
%   Dimensiona in automatico la sezione in base ai parametri e limiti
%   predefiniti, in funzione della sollecitazione agente. Sono previsti al
%   massimo 2 righe di barre allo stesso livello di diametro differente.
%   Una popolazione di risultati viene generata in funzioni delle variabili
%   diametro delle barre, numero di barre per riga, numero di righe i cui
%   limiti sono definiti rispettivamente in fi_lim, Nb1 e Nb2. Questi
%   3 parametri sono matrici Nx2 dove N = lentgh(d);
%   Il dimensionamento è effettuato in condizioni di armatura semplice,
%   i.e. è ignorato il contributo dell'armatura in compressione.
%   Variabili di input:
%       sezione:    matrice Nx4, dove la N riga è composta dai 4 elementi
%       che caratterizzano il rettangolo (xm, ym, dx, dy);
%       armatura: [d, nb, fi] matrice Nx3 che contiene i dati base della
%           armatura. d:  altezza utile della sezione dal lembo compresso (il lembo compresso coincide con l'asse delle ascisse)
%           nb: numero di barre per livello fi: diametro delle armature di base
%       fi_lim: [fi_min, fi_max] coppia di valori che indicano il limite
%       inferiore ed il limite superiore dei diametri da utilizzare per
%       l'analisi;
%       Nb1:    vettore 1xN contenente il numero di barre possibili
%       Nb2:    come Nb1, ma per il diametro 2 (porre Nb2 = 0 per nessuna
%       barra aggiuntiva, oppure 1:x per barre aggiuntive da 0 a x, perché
%       la possibilità Nb2 = 0 è sempre presa in considerazione.
%       fiv:    vettore 1xN contenente i diametri possibili
%       fck:    valore di resistenza caratteristica del cls
%       Med:    momento sollecitante per il dimensionamento (>0)
%       Ned:    sforzo normale per il dimensionamento (>0 significa
%       compressione)
%   Argomenti opzionali:
%       tipo: 'elastica', 'plastica'; argomento da passare alla funzione
%       calcoloNM, ne determina il tipo di analisi da effettuare
%       lock: 'no', 'yes'; se lock='yes' allora Nb2(i)=Nb1(j) in ogni caso, non
%       ci sono combinazioni differenti (pertanto Nb2 diventa obsoleto come
%       valore di input). Se lock='no' l'analisi comprende le diverse
%       combinazioni di barre per cui Nb2(i)~=Nb1(j);
%       precisione: 12; specifica la precisionde del calcolo del Mrd

%% Valori di default ed estrazione argomenti opzionali
tipo = 'elastica';
lock = false;
precisione = 12;
ottimizzazione = true;

loadOptionalArgument(varargin);

%% funzioni e valori di default
As.fun = @(n,fi) n*pi*fi^2/4; % area di armatura per il diametro e numero di barre
func = ['calcoloNM(x, sezione, d, As.tot, def_not, mat.cls.f_cd, mat.steel.f_yd, ''' tipo ''');'];
x = 0;
d = armatura(:,1);  

%% estrazione dei vettori della matrice "sezione"
[~, ~, ~, ~, B, H, uym] = estrazioneSezione(sezione);

%% estrazione delle deformazioni notevoli
mat.cls = derivaCaratteristicheCA(fck);
mat.steel = derivaCaratteristicheAcciaio;
defstrcls = {'ecu','ec2','ec3'};
defstrsteel = {'esu','eyd'};
for i = 1:length(defstrcls)
    def_not.(defstrcls{i}) = mat.cls.(defstrcls{i});
end
for i = 1:length(defstrsteel)
    def_not.(defstrsteel{i}) = mat.steel.(defstrsteel{i});
end


%% armatura minima
d_lim = max(d);   %
a = find(uym <= d_lim, 1, 'last');   % estremo inferiore dell'intorno di d
b = find(uym >= d_lim, 1, 'first');  % estremo superiore dell'intorno di d
B_medio = (B(a,1) + B(b,1))/2;  % larghezza media della sezione all'altezza dell'armatura
As.min = 0.26 * mat.cls.f_ctm/mat.steel.f_yk * B_medio * d_lim;  % minimo di armatura secondo normativa, per l'estremo superiore ed inferiore

%% determinazione barre minime
if lock
    length_Nb2 = 1;
    Nb2 = Nb1;
else
    length_Nb2 = length(Nb2);
end

%% determina se deve essere verificato il momento negativo oppure positivo
if Med < 0
    % se il momento è negativo, si inverte l'asse delle ordinate, in quanto
    % l'origine coincide con il lembo compresso della sezione.
    sezione(:,2) = H - sezione(:,2);    % inversione delle ym
    armatura(:,1) = H - armatura(:,1) ; % inversione delle altezze utili d
end

% estrae d e riordina gli elementi dal più piccolo al più grande.
% L'elemento più grande corrisponde alla riga che si va a dimensionare. Con
% idex_d si riordina la matrice armatura in modo tale da avere
% corrispondenza tra altezze utili d e barre di armatura.
[d, index_d] = sort(armatura(:,1)); 
armatura = armatura(index_d,:);

% calcolo dell'area di armatura invariante
As.tot = zeros(size(d));
for i_d = 1:length(d)-1
    As.tot(i_d) = As.fun(armatura(i_d, 2), armatura(i_d, 3));  % nb x fi
end
arm = struct(...
    'd', armatura(:,1),...
    'nb1', armatura(index_d,2),...
    'nb2', zeros(size(d)),...
    'fi1', armatura(index_d,3),...
    'fi2', zeros(size(d)),...
    'As1', As.tot(index_d),...
    'As2', zeros(size(d))...
    );

%% Inizializzazione tabella dei risultati
ris = struct('nb1',0, 'nb2',0, 'fi1', 0, 'fi2', 0, 'As1', 0, 'As2', 0, 'As_tot', 0, 'Mrd_fi1', 0, 'Mrd', 0, 'ratio', 0, 'armatura', arm);
ris = repmat(ris, length(Nb1)*length_Nb2*length(fiv), 1);
i_r = 0;

%% Impostazioni per ottimizzazione == false
% Nel caso non si voglia eseguire il ciclo di ottimizzazione, si riduce il
% vettore fiv ad 1 elemento in modo tale che coincida con il diametro per
% le barre fi1. Il secondo elemento di fiv è invece il diametro per fi2.
if ~ottimizzazione
    fi2 = fiv(2);
    fiv = fiv(1);
end

%% inizio ciclo principale
for i_nb1 = 1:length(Nb1)
    fi = zeros(2,1);    % vettore 2x1: la componente 1 è relativo al diametro 1, la componente 2 relativa al diametro 2
    
    for i_fiv = 1:length(fiv)
        % reset di fi
        fi(1) = fiv(i_fiv); % itero sui diametri base
        fi(2) = 0;  % azzero il secondo diametro ### (non dovrebbe servire)
        
        for i_nb2 = 1:length_Nb2
            if lock
                i_nb2locked = i_nb1;
            else
                i_nb2locked = i_nb2;
            end
            
            Mrd = 0;
            j = 0;
            % il ciclo determina il diametro 2 in combinazione con il
            % diametro 1 per cui Mrd > Med. Il diametro 1 non varia nel
            % ciclo, solo il 2.
            while Mrd <= abs(Med)
                
                if ottimizzazione
                    % aggiornamento del diametro 2
                    fi(2) = sign(j) * (fiv(i_fiv) + 2 *(j-1));
                    
                    % condizione di interruzione ciclo
                    if fi(2) > max(fiv)
                        break
                    end
                else
                    fi(2) = j*fi2;
                end
                
                % calcolo dell'area di armatura
                As.f1 = As.fun(Nb1(i_nb1), fi(1));
                As.f2 = As.fun(Nb2(i_nb2locked), fi(2));
                As.tot(end) = As.f1 + As.f2;
                if As.f2 > 0 && As.tot(end) < As.min
                    Mrd = 0;
                else
                    % calcolo del momento resistente
                    brent('x','[N, M]', func, Ned, 0, H, 'precisione', precisione);
                    Mrd = M*1e-6;
                end
                
                if fi(2) == 0
                    Mrd_fi1 = Mrd;
                end
                
                % se nel caso in cui non ci siano barre aggiuntive, allora
                % termina il ciclo perché non ha soluzione
                if Nb2(i_nb2locked) == 0
                    break
                end
                
                % condizione di fine ciclo per caso specifico
                if ~ottimizzazione
                    % se fi2 > 0, esegue sempre 2 cicli. Atrimenti
                    % interrompe il ciclo alla seconda iterazione.
                    if j == 0 && fi2 ~= 0
                        Mrd = 0;    
                    else
                        % se non si esegue il ciclo di ottimizzazione, si esce
                        % sempre alla prima iterazione
                        break
                    end
                end
                
                % aggiornamento variabile di ciclo
                j = j + 1;
            end
            
            % salvataggio dei risultati
            if Mrd > abs(Med)
                i_r = i_r + 1;
                ris(i_r).nb1 = Nb1(i_nb1);
                ris(i_r).fi1 = fi(1);
                if fi(2) ~= 0
                    ris(i_r).nb2 = Nb2(i_nb2locked); % Se il diametro è nullo è inutile riportare una quantità di barre pertanto è lasciato pari a 0
                end
                ris(i_r).fi2 = fi(2);
                ris(i_r).As1 = As.f1;
                ris(i_r).As2 = As.f2;
                ris(i_r).As_tot = As.tot(end);
                ris(i_r).Mrd_fi1 = Mrd_fi1;
                ris(i_r).Mrd = Mrd;
                ris(i_r).ratio = abs(Med)/Mrd;
                % aggiorno la variabile arm con il valore corrente
                if Med < 0
                    i_obb = 1;
                else
                    i_obb = length(d);
                end
                arm.nb1(i_obb) = ris(i_r).nb1;
                arm.nb2(i_obb) = ris(i_r).nb2;
                arm.fi1(i_obb) = ris(i_r).fi1;
                arm.fi2(i_obb) = ris(i_r).fi2;
                arm.As1(i_obb) = ris(i_r).As1;
                arm.As2(i_obb) = ris(i_r).As2;
                ris(i_r).armatura = arm;
            end
            
            % condizione di fine ciclo
            if and(or(fi(2) == 0, fi(1) == fi(2)), Mrd > abs(Med))
                break
            end
            
        end
        if and(fi(2) == 0, Mrd > abs(Med))
            break
        end
    end
end

%% eliminazione dei campi vuoti
while ~isempty(ris) && ris(end).nb1==0
    ris(end) = [];
end

end