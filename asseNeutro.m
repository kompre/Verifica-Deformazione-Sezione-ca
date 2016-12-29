function [ x, A, S, J ] = asseNeutro( b, h, As, d, N, M ,n)
%ASSENETURO calcola asse neutro sezione rettangolare soggetta a flessione o
%pressoflessione
%   calcola la posizione dell'asse neutro rispetto al lembo compresso di
%   cls per sezione rettangolare soggetta a flessione o pressoflessione
%   retta
%   b: base sezione
%   h: altezza sezione
%   As: area di acciaio
%   d: distanza delle barre di acciaio dal lembo compresso
%       As e d devono avere la stessa dimensione
%   N: sforzo normale
%   M: momento flettente
%   n: rapporto tra modulo elastico acciaio/cls

if N == 0
    x0 = h/2;    % variabile di primo tentativo
    err = 1;    % inzializzazione variabile di controllo del ciclo
    while err > 1e-6   %precisione desiderata
        S = b*x0^2/2 + n*sum(As.*d);    % momento statico rispetto al lembo compresso di cls;
        A = b*x0 + n*sum(As);  % area omogeneizzata
        x = S/A;    % posizione dell'asse neutro al passo i+1;
        err = abs(x-x0); % errore tra due iterazioni successive
        x0 = x; %aggiornamento della variabile di iterazione
    end
else
    ecc = M/N;  %eccentricità del punto di applicazione
    u = ecc-h/2;    % distanza del punto di applicazione dal lembo compresso
    x0 = h/2;
    err = 1;
    while err > 1e-6
        yp = x0+u;  % distanza dell'asse neutro dal punto di applicazione
        S = 1/2*b*(yp^2-u^2) + n*sum(As.*(d+u)); % momento statico rispetto al punto di applicazione
        I = 1/3*b*(yp^3-u^3) + n*sum(As.*(d+u).^2); %momento di inerzia rispetto al punto di applicazione
        yp1 = I/S; % distanza dell'asse neutro dal punto di applicazione dell'iterazione successiva
        x = yp1-u; % posizione asse neutro rispetto al lembo compresso
        err = abs(x-x0);    %errore tra due iterazioni successive
        x0 = x; % aggiornamento variabile di iterazione
    end
end
A = b*x + n*sum(As);  % area omogeneizzata
S = b*x^2/2 + n*sum(As.*(x-d)); % momento statico omogeneizzato rispetto all'asse neutro
J = b*x^3/3 + n*sum(As.*(d-x).^2);  % momento di inerzia omogeneizzato rispetto all'asse neutro
end

