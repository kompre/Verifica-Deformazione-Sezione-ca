function [ prop ] = rebarProperty( armatura )
%SECTIONPROPERTY calcolo delle proprietà inerziali per barre di acciaio
%   calcolo delle proprietà inerziali per le barre di armatura. Le barre di
%   armatura si considerano come area puntiforme, pertanto il calcolo del
%   momento di inerzia è semplificato in A*(d-yg)^2
%       1. xm: coordinata del baricentro della barra
%       2. ym: coordinata del baricentro della barra
%       3. As: Area dell'armatura

%% estrazione dati sezione
xm = armatura(:,1);  % coordinata del baricentro del rettangolo i-esimo
ym = armatura(:,2);  % coordinata del baricentro del rettangolo i-esimo
As = armatura(:,3);  % larghezza del rettangolo i-esimo


%% Calcolo proprietà
prop.A = sum(As);    % somma delle aree di armatura
prop.Sx = sum(As.*ym);  % momento statico rispetto all'asse x (y = 0)
prop.Sy = sum(As.*xm);  % momento statico rispetto all'asse y (x = 0)
%% coordinate baricentro sezione
prop.yg = prop.Sx/prop.A; 
prop.xg = prop.Sy/prop.A;  
%% momento di inerzia rispetto al baricentro della sezione
prop.Jxg = sum(As.*(ym - prop.yg).^2);   % momento di inerzia rispetto all'asse x passante per il baricentro (y = yg)
prop.Jyg = sum(As.*(xm - prop.xg).^2);   % momento di inerzia rispetto all'asse y passente per il baricentro (x = xg)



end

