function loadOptionalArgument( argument )
%LOADOPTIONALARGUMENT carica nella funzione chimante le opzioni specificate
%in varargin in coppia 'key', 'value'. entrambi gli argomenti devono essere
%delle stringhe
%   Detailed explanation goes here

while ~isempty(argument)
    evalin('caller', [argument{1} ' = ' argument{2} ';']);
    argument(1:2) = [];
end

end

