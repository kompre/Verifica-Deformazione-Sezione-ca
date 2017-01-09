function [ rebar ] = rebarDistr( nb, diam, d, a, b )
%REBARDIST distribuisce le barre di armatura sull'intervallo di estremi
% [a,b]
%   Detailed explanation goes here



rebar = zeros(sum(nb),3);
k = 0;
for i = 1:length(d)
    dxs = (b-a)/nb(i); % interasse tra due armature
    xs = a+dxs/2:dxs:b; % vettore delle ascisse delle barre
    for j = 1:length(xs)
        k = k+1;
        rebar(k,:) = [xs(j), d(i), pi*diam(i)^2/4];
    end
end

end