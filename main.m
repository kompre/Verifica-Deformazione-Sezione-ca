clearvars
clc
%% dati sezione
% dati geometrici
b = [1000]
h = [300]
% dati armatura
As = [565;565];
d = [40;260];
% dati materiale
E.cls = 31e3;
E.steel = 210e3;
n = E.steel/E.cls; 
% dati sollecitazione
N = 0;
M = 51.8;
%% calcolo carettristiche sezione c.a. (rispetto al baricento della sezione)
A.cls = b*h;
S.cls = 0;
J.cls = b*h^3/12;
%% calcolo inerzia sezione interamente reagente
sez = rettangolo(b,h,[0,0],[0,300],100,100);
Ax = (sez(:,3)'*sez(:,4))
Sx = (sez(:,3).*sez(:,4))'*sez(:,2)
Jx = sum(sez(:,3).*sez(:,4).*(sez(:,4).^2/12 + sez(:,2).^2))
%%
yg = Sx/Ax
J0 = Jx - Ax*yg^2




%% calcolo inerzia sezione parzializzata
[x, A.omog, S.omog, J.omog] = asseNeutro(b, h, As, d, N, M, n);
Jg = J.omog + A.omog * (h/2-x)^2;
alpha.fess = (E.cls*J.cls) / (E.steel/n*Jg)