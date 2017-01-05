function [ ris ] = dimenSezione( sezione, d, fi_lim, Nb1, Nb2, mat, Med, Ned, varargin )
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

%% Estrazione argomenti opzionali
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'tipo'
            tipo = varargin{2};
        case 'lock'
            lock = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

if ~exist('tipo', 'var')
    tipo = 'elastica';
end

if ~exist('lock', 'var')
    lock = 'no';
end


%% funzioni
As.fun = @(n,fi) n*pi*fi^2/4; % area di armatura per il diametro e numero di barre
func = ['calcoloNM(x,sezione,d,As.tot,def_not,mat.cls.f_cd,mat.steel.f_yd,''' tipo ''');'];
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

x = 0;

%% estrazione delle deformazioni notevoli
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
fiv = fi_lim(1):2:fi_lim(2); % vettore dei diametri di armatura possibili
nb1v = Nb1(1):Nb1(2);  % vettore del numero di barre possibili

if strcmp(lock,'no')
    nb2v = Nb2(1):Nb2(2);  % vettore del numero di righe possibili
    nb2v_ll = nb2v;
    nb2v_ll = 1;
    nb2v = nb1v;
end

% Inizializzazione tabella dei risultati
ris = struct('nb1',0, 'fi1', 0,'nb2',0, 'fi2', 0, 'As1', 0, 'As2', 0, 'As_tot', 0, 'Mrd', 0, 'ratio', 0);
ris = repmat(ris, length(nb1v)*length(nb2v)*length(fiv), 1);
i_r = 0;

for i_nb1 = 1:length(nb1v)
    fi = zeros(2,1);    % vettore 2x1: la componente 1 è relativo al diametro 1, la componente 2 relativa al diametro 2
    
    for i_fiv = 1:length(fiv)
        % reset di fi
        fi(1) = fiv(i_fiv); % itero sui diametri base
        fi(2) = 0;  % azzero il secondo diametro ### (non dovrebbe servire)
        
        
        for i_nb2 = 1:length(nb2v_ll)
            if strcmp(lock,'yes')
                i_nb2locked = i_nb1;
            else
                i_nb2locked = i_nb2;
            end
            
            Mrd = 0;
            j = 0;
            % il ciclo determina il diametro 2 in combinazione con il
            % diametro 1 per cui Mrd > Med. Il diametro 1 non varia nel
            % ciclo, solo il 2.
            while Mrd < abs(Med)
                
                % aggiornamento del diametro 2
                fi(2) = sign(j) * (fiv(i_fiv) + 2 *(j-1));
                
                % condizione di interruzione ciclo
                if sum(fi(2) > max(fi_lim))
                    break
                end
                
                % calcolo dell'area di armatura
                As.f1 = As.fun(nb1v(i_nb1),fi(1));
                As.f2 = As.fun(nb2v(i_nb2locked),fi(2));
                As.tot = As.f1 + As.f2;
                if As.tot < As.min
                    Mrd = 0;
                else
                    % calcolo del momento resistente
                    dicotomico('x','[N, M]', 'calcoloNM(x,sezione,d,As.tot,def_not,mat.cls.f_cd,mat.steel.f_yd,''elastica'');', Ned, 0, H, 6)
                    Mrd = M*1e-6;
                end
                
                % se nel caso in cui non ci siano barre aggiuntive, allora
                % termina il ciclo perché non ha soluzione
                if nb2v(i_nb2locked) == 0
                    break
                end
                
                % aggiornamento variabile di ciclo
                j = j + 1;
            end
            %% salvataggio dei risultati
            if Mrd >= Med
                i_r = i_r + 1;
                ris(i_r).nb1 = nb1v(i_nb1);
                ris(i_r).fi1 = fi(1);
                if fi(2) ~= 0
                    ris(i_r).nb2 = nb2v(i_nb2locked); % Se il diametro è nullo è inutile riportare una quantità di barre pertanto è lasciato pari a 0
                end
                ris(i_r).fi2 = fi(2);
                ris(i_r).As1 = As.f1;
                ris(i_r).As2 = As.f2;
                ris(i_r).As_tot = As.tot;
                ris(i_r).Mrd = Mrd;
                ris(i_r).ratio = abs(Med)/Mrd;
            end
            %% condizione di fine ciclo
            if or(fi(2) == 0, fi(1) == fi(2))
                break
            end
            
        end
        if fi(end) == 0
            break
        end
    end
end
%% eliminazione dei campi vuoti
risTable = struct2table(ris);
campi_pieni = risTable.fi1>0;    % vettore logico dei campi non vuoti
risTable = risTable(campi_pieni, :);
ris = table2struct(risTable);   % riconverto la tabella in struct

end