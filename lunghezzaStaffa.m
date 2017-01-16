function [ Lst, peso ] = lunghezzaStaffa( H, B, nb, copriferro, varargin)
%LUNGHEZZASTAFFA calcolo della lunghezza della singola staffa, data la
%sezione ed il copriferro
%   Determina la lunghezza della staffa rettangolare per una sezione di
%   dimensioni BxH, in funzione del copriferro e del numero di bracci nb.
%   Si assume che la staffa abbia forma rettangolare a due bracci, pertanto
%   il nb determina il numero di staffe rettangolari di calcolo.
%   Inoltre si assume che i bracci delle staffe siano equidistantti tra
%   loro.

%% load optional argument
peso_acciaio = 7850;
fi = 0;
loadOptionalArgument(varargin)
%% estrazione del copriferro

lato = {'inf', 'sup', 'sx', 'dx'};
for i = 1:4
    if length(copriferro) > 1
        c.(lato{i}) = copriferro(i);
    else
        c.(lato{i}) = copriferro;
    end
end

%% determinazione del numero di staffe rettangolari
nSt = nb/2;

%% dimensioni nette delle staffe
b.netto = B - (c.sx + c.dx);
h.netto = H - (c.inf + c.sup);
b.tot = b.netto*(1+(nSt-1)/(nb-1));  % dimensione cumulata delle staffe che tiene conto della sovrapposizione

%% Calcolo lunghezza della staffa
Lst = nb*h.netto + 2*b.tot;

%% Calcolo del peso della staffa
if fi > 0
    Area = pi*fi^2/4;  % area dell'armatura
    peso = Area*1e-6 * Lst*1e-3 * peso_acciaio;
end


end

