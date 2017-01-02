clearvars
clc
%% dati sezione
% dati geometrici
b = [1000];
h = [300];
l = 4;  % [m] luce libera di inflessione
% dati armatura
nb = [5;5];
diam = [12;12];
d = [40;260];
% dati materiale
mat.cls = derivaCaratteristicheCA(25,30);
mat.steel = derivaCaratteristicheAcciaio;
n = mat.steel.Es/mat.cls.E_cm; % coefficiente di omogeneizzazione
% dati sollecitazione
Fd = 66.3;  % [kN/m2] carico agente in condizioni SLE
Ned = 0;
M = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;   % funzione del momento sollecitante per trave appoggio appoggio
% dati punti di calcolo
dx = .001;  % passo dell'integrazione
x = 0:dx:l; % vettore delle ascisse dei punti di calcolo
%% calcolo inerzia sezione interamente reagente
sez = rettangolo(b,h,0,0,100,100);    % discretizzazione della sezione in c.a. in rettangoli infitesimali
reb = rebarDistr(nb, diam, d, -b/2, b/2);   % distribuzione delle barre di armatura nella sezione

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
W_min = prop.ca.Jxg/min(prop.ca.yg,(h-prop.ca.yg));
Mf = W_min*mat.cls.f_ctm*1e-6; % momento di prima fessurazione [kNm]
alpha.nf = (mat.cls.E_cm * prop.cls.Jxg) / (mat.steel.Es/n *prop.ca.Jxg); % rapporto momento di inerzia sezione non fessurata
%% Calcolo asse neutro per sezione fessurata
[sezfes.x, sezfes.A, sezfes.Sx, sezfes.Jx] = asseNeutro(b,h,reb(:,3),reb(:,2),Ned,Mf,n); % N e M sono necessari sono nel caso di N > 0

Med = M(Fd,l,x);
curvatura = zeros(size(Med));
% for i = 1:length(Med)
%     if abs(Med(i)) <= Mf
%         J = prop.ca.Jxg;
%     else
%         J = sezfes.Jx;
%     end
%     curvatura(i) = Med(i)*1e6/(mat.cls.E_cm*J);
% end
curvatura = Med*1e6/(mat.cls.E_cm*prop.ca.Jxg);

rotazione = zeros(size(curvatura));
spostamento = zeros(size(curvatura));

rotazione = integraleIndefinito(curvatura,x*1e3,0);
spostamento = integraleIndefinito(rotazione,x*1e3,0);

%%
    


figure(1)
plot(x,Med)
figure(2)
hold on
plot(x,curvatura)
plot(x,rotazione)
plot(x,spostamento)

