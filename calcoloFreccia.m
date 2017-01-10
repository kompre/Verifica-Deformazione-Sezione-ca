function [s_min] = calcoloFreccia( L, geom, arm, fck, c, dx, soll, tipo, valueTipo)
%  Calcolo della freccia elastica della trave in ca
%   calcolo della freccia elastica per trave in c.a. di luce 'L', di
%   sezione 'geom' e armatura 'arm', soggetta alla sollecitazione 'soll'.
%   La funzione offre 2 tipologie di calcolo, specificate in varargin:
%       'accurato': calcola la deformata tramite integrale indefinito a
%       partire dalla curvatura. valueTipo = [rot, spo] i.e. costanti
%       dell'integrale per rotazione (1) e spostamento (2);
%       'semplificato': calcola la freccia con la funzione 'f_el'
%       specificata in valueTipo;

%% estrazione dati
cost.c = c;     % costante per il calcolo del parametro zeta (0.5, 1)

% tipologia di analisi
% se tipo='accurato', allora valueTipo è un vettore 2x1 dove il primo
% elemento è la costante di integrazione della rotazione ed il secondo
% elemento è la costante di intergazione dello spostamento;
% se tipo='semplificato', allora valueTipo è una funzione anonima, che
% specifica il valore della freccia max in mezzeria della trave;
switch lower(tipo)
    case 'accurato'
        cost.rot = valueTipo(1); 
        cost.spo = valueTipo(2);
    case 'semplificato'
        tipo = 'semplificato';
        f_el = valueTipo;
end

% dati sollecitazione
Fd = soll.q;  % [kN/m2] carico agente in condizioni SLE
Ned = soll.N;
M = soll.M;   % funzione del momento sollecitante per trave 1 campata

% dati vari
x = 0:dx:L; % vettore delle ascisse dei punti di calcolo
mat.cls = derivaCaratteristicheCA(fck);
mat.steel = derivaCaratteristicheAcciaio;
n = mat.steel.E_s/mat.cls.E_cm; % coefficiente di omogeneizzazione

%% Proprietà della sezione
sez = rettangolo(geom.b, geom.h, geom.x0, geom.y0, 100, 100);    % discretizzazione della sezione in c.a. in rettangoli infitesimali
reb = rebarDistr(arm.nb, arm.diam, arm.d, -geom.b/2, geom.b/2);   % distribuzione delle barre di armatura nella sezione

% caratteristiche sezione cls
prop.cls = sectionProperty(sez);

% caratteristiche sezione acciaio
prop.steel = rebarProperty(reb);

% caratteristiche sezione ca (cls+acciaio)
prop.ca.A = prop.cls.A + n*prop.steel.A;    % area omogeneizzata rispetto al cls
prop.ca.Sx = prop.cls.Sx + n*prop.steel.Sx; % momento statico omogeneizzato rispetto al cls
prop.ca.yg = prop.ca.Sx/prop.ca.A;  % baricentro omogeneizzato sezione c.a.
prop.ca.Jxg = prop.cls.Jxg + prop.cls.A*(prop.cls.yg-prop.ca.yg)^2 + n*(prop.steel.Jxg + prop.steel.A*(prop.steel.yg-prop.ca.yg)^2);

%% Calcolo del momento di prima fessurazione
W_min = prop.ca.Jxg/min(prop.ca.yg,(sum(geom.h)-prop.ca.yg));
Mf = W_min*mat.cls.f_ctm*1e-6; % momento di prima fessurazione [kNm]
alpha.nf = (mat.cls.E_cm * prop.cls.Jxg) / (mat.steel.E_s/n *prop.ca.Jxg); % rapporto momento di inerzia sezione non fessurata
%% Calcolo asse neutro per sezione fessurata
[sezfes.x, sezfes.A, sezfes.Sx, sezfes.Jx] = asseNeutro(geom.b, geom.h, reb(:,3),reb(:,2),Ned,Mf,n); % N e M sono necessari sono nel caso di N > 0

%% Calcolo delle deformazioni
% La deformata è stimata moltiplicando la freccia elastica per un fattore
% amplificativo pari al rapporto tra il momento di inerzia della sezione
% lorda di cls (è ignorato il contributo delle barre) e il momento di
% inerzia intermedio tra la sezione non fessurata e la sezione fessurata
% (con barre).

Med = M(Fd,L,x);    %calcolo del momento sollecitante lungo l'asse

momenti_notevoli = [Med(1), max(Med), Med(end)];    % momenti alle sezioni notevoli
j = 0;
for i = 1:length(momenti_notevoli)
    beta = Mf/momenti_notevoli(i);
    if beta < 1
        j = j + 1;
        zeta(j) = 1 - cost.c*beta^2;
    end
end
Jef = zeta * sezfes.Jx + (1 - zeta) * prop.cls.Jxg; % momento di inerzia efficace
Jef = mean(Jef);

if strcmp(tipo, 'accurato')
    curvatura = Med*1e6/(mat.cls.E_cm * Jef);  % curvatura
    rotazione = integraleIndefinito(curvatura, x*1e3, cost.rot);  % rotazione
    spostamento = integraleIndefinito(rotazione, x*1e3, cost.spo);  % rotazione
    s_min = min(spostamento); % freccia massima
elseif strcmp(tipo, 'semplificato')
    alpha.fes = prop.cls.Jxg/Jef;   % parametro di amplificazione della freccia elastica
    s_min = alpha.fes * f_el(Fd, L*1e3, mat.cls.E_cm, prop.cls.Jxg);
end

end

