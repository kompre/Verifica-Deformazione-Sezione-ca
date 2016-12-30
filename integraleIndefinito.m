function [ Fx ] = integraleIndefinito( fx, x, cost )
% Calcolo dell'intergrale indefinito di fx, nell'intervallo x, a meno di
% cost
%   Il valore dell'integrale indefinito è ottenuto sommando i valori delle
%   aree infinitesimi nell'intervallo [0,i], dove i è la componente i-esima
%   del vettore fx. Ne risulta un vettore Fx che ha la stessa dimensione di
%   fx. Poiché il metodo di calcolo è definito sull'intervallo [0,i] è
%   necessario calcolare la costante c0 per traslare il vettore delle
%   ordinate della quantità necessaria affinché le ascisse combacino con le
%   ordinate.
%   Il metodo di integrazione utilizzato è quello dei trapezi, pertanto di
%   misura il valore medio di fx tra due punti consecutivi. Tale valore
%   (fm) è moltiplicato per il passo di integrazione dx e costituisce
%   l'area infinitesimale al passo i-esimo. 
%       fx: un vettore contenente i valori della funzione integranda;
%       x: vettore delle ascisse (valori ordinati)
%       cost:   costante di integrazione

dx = x(2) - x(1);     % passo di integrazione ( si assume che x sia equispaziato)
a = min(x); % limite inferiore dell'intervallo x
xa = min(0,a):dx:max(0,a);
dim(1) = length(x);
dim(2) = length(xa);
f = struct;
for j = 1:2
    f(j).fm = zeros(1,dim(j)-1); % NB: fm ha dimensione N-1, dove N è la lunghezze di fx
    f(j).fmdx = zeros(size(f(j).fm));    % fmdx è un vettore di dimensioni uguali a fm
    for i = 2:dim(j)
        f(j).fm(i) = (fx(i) + fx(i-1))/2;    % valore intermedio della funzione nell'intervallo (metodo dei trapezi)
        f(j).fmdx(i) = f(j).fm(i)*dx;
    end
end

Fx = zeros(size(fx));
if length(xa)~=1
    c0 = sum(f(2).fmdx(1:dim(2)));
else
    c0 = 0;
end

Fx(1) = sign(a)*c0 + cost;   % inizializzo il primo valore dell'integrale
for i = 1:length(fx)-1
    Fx(i+1) = sum(f(1).fmdx(1:i)) +sign(a)*c0 + cost;    
end
% Fx(1) = cost;   % inizializzo il primo valore dell'integrale
% for i = 1:length(fx)-1
%     Fx(i+1) = sum(f(1).fmdx(1:i)) + cost;    
% end

end

