function [ xm, ym, db, dh, B, H, uym] = estrazioneSezione( sezione )
%ESTRAZIONESEZIONE Summary of this function goes here
%   Detailed explanation goes here

xm = sezione(:,1);  % il primo campo alla coordinata in xm del baricentro dei rettangoli
ym = sezione(:,2);  % il secondo campo è riservato alla coordinata in ym del baricentro dei rattngoli
db = sezione(:,3);  % il terzo campo è riservato alla larghezza in x dei rettangoli
dh = sezione(:,4);  % il quarto campo è riservato alla larghezza in y dei rettangoli

[uym, ih, ~] = unique(ym);    % elimina le componenti duplicate
H = sum(dh(ih));    % massima altezza della sezione (somma solo una volta dh alla quota ym)
B = zeros(length(uym),2);
for i = 1:length(uym)
    B(i,:) = [sum( abs( db( ym == uym(i) ) ) ), uym(i)]; % larghezza della sezione in funzione dell'altezza
end



end

