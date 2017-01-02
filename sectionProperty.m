function [ prop ] = sectionProperty( sezione )
%SECTIONPROPERTY calcolo delle proprietà principali di una sezione
%   calcolo delle proprietà fondamentali di una sezione qualunque.
%   Sezione è una matrice Nx4, dove N rappresenta il numero di rettangoli
%   i-esimi in cui è discretizzara la sezione. Per ciascun rettangolo si
%   specifica:
%       1. xm: coordinata del baricentro del rettangolo i-esimo rispetto
%       all'origine degli assi
%       2. ym: coordinata del baricentro del rettangolo i-esimo rispetto
%       all'oridine degli assi
%       3. dx: dimensione del rettangolo in direzione x
%       4. dy: dimensione del rettangolo in direzione y

%% estrazione dati sezione
xm = sezione(:,1);  % coordinata del baricentro del rettangolo i-esimo
ym = sezione(:,2);  % coordinata del baricentro del rettangolo i-esimo
dx = sezione(:,3);  % larghezza del rettangolo i-esimo
dy = sezione(:,4);  % altezza del rettangolo i-esimo

%% Calcolo proprietà
prop.A = dx'*dy;    % somma delle aree (1xN)x(Nx1)=(1x1)
prop.Sx = sum(dx.*dy.*ym);  % momento statico rispetto all'asse x (y = 0)
prop.Sy = sum(dx.*dy.*xm);  % momento statico rispetto all'asse y (x = 0)
%% coordinate baricentro sezione
prop.yg = prop.Sx/prop.A; 
prop.xg = prop.Sy/prop.A;  
%% momento di inerzia rispetto al baricentro della sezione
prop.Jxg = sum(dx.*dy.^3/12 + dx.*dy.*(ym - prop.yg).^2);   % momento di inerzia rispetto all'asse x passante per il baricentro (y = yg)
prop.Jyg = sum(dy.*dx.^3/12 + dx.*dy.*(xm - prop.xg).^2);   % momento di inerzia rispetto all'asse y passente per il baricentro (x = xg)



end

