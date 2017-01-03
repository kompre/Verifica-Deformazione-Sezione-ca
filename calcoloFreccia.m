function [s_min] = calcoloFreccia( L, geom, arm, mat, soll, cost, dx)
%% dati sezione
% dati geometrici
% b = [1000];
% L = 4;  % [m] luce libera di inflessione
% dati armatura
% nb = [5;5];
% diam = [12;12];
% d = [40;260];
% dati materiale
% mat.cls = derivaCaratteristicheCA(25,30);
% mat.steel = derivaCaratteristicheAcciaio;
n = mat.steel.Es/mat.cls.E_cm; % coefficiente di omogeneizzazione
% dati sollecitazione
Fd = soll.q;  % [kN/m2] carico agente in condizioni SLE
Ned = soll.N;
M = soll.M;   % funzione del momento sollecitante per trave 1 campata
% cost.c = 0.5;    % parametro per il calcolo di zeta [0.5 per carichi lunga durata/ripetuti, 1 per breve duarat]
% % condizioni al contorno
% cost.rot = 0;   % c.i. rotazione
% cost.spo = 0;   % c.i. spostamento
% dati punti di calcolo
% dx = .001;  % passo dell'integrazione
x = 0:dx:L; % vettore delle ascisse dei punti di calcolo
%%
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
alpha.nf = (mat.cls.E_cm * prop.cls.Jxg) / (mat.steel.Es/n *prop.ca.Jxg); % rapporto momento di inerzia sezione non fessurata
%% Calcolo asse neutro per sezione fessurata
[sezfes.x, sezfes.A, sezfes.Sx, sezfes.Jx] = asseNeutro(geom.b, geom.h, reb(:,3),reb(:,2),Ned,Mf,n); % N e M sono necessari sono nel caso di N > 0

%% Calcolo delle deformazioni
% La deformata � stimata moltiplicando la freccia elastica per un fattore
% amplificativo pari al rapporto tra il momento di inerzia della sezione
% lorda di cls (� ignorato il contributo delle barre) e il momento di
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

curvatura = Med*1e6/(mat.cls.E_cm * Jef);  % curvatura
rotazione = integraleIndefinito(curvatura, x*1e3, cost.rot);  % rotazione
spostamento = integraleIndefinito(rotazione, x*1e3, cost.spo);  % rotazione
s_min = min(spostamento); % freccia massima
end

