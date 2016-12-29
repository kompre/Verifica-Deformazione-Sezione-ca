function [ Fx ] = integraleIndefinito( fx, x, cost )
% Calcolo dell'intergrale indefinito della funzione fx, nella variabile var
% Calcola l'integrale indefinito per ogni elemento del vettore fx, a meno della costante cost (default: cost = 0)
% x è il vettore delle ascisse corrispondente e serve per determinare il dx
% Fx è vettore dei valori integrali

fm = zeros(1,length(fx)-1); % NB: fm ha dimensione N-1, dove N è la lunghezze di fx
fmdx =zeros(size(fm));
for i = 2:length(x)
    dx = x(i) - x(i-1);     % ampiezza dell'intervallo tra due punti successivi
    fm = fx(i) - fx(i-1);    % valore intermedio della funzione nell'intervallo (metodo dei trapezi)
    fmdx = fm*dx;   % 
end

Fx = zeros(size(fx));
Fx(1) = cost;   % inizializzo il primo valore dell'integrale
for i = 1:length(fx)-1
    Fx(i+1) = sum(fmdx(1:i)) + cost;
end
plot(x,Fx)
Fx(end)

