% Caricamento dati di default
close all
pre
save('data_pre.mat')
var2save = {'arm', 'cost', 'cu', 'f_lim', 'geom', 'L', 'num', 'peso_acciaio', 'RZ', 'cedimento_globale', 'sollecitazioni', 'k_Minf', 'k_Msup'};
var2var = {...
    %non compattabili
    'sollecitazioni.pl.q', {'50 - 5.6'};...
    'sollecitazioni.pl.qSLU', {'1.3 * 8.8 + 1.5 * 50 - 24'};...
    'num.campate' , {'16'};...
    'num.pali_tr', {'15'};...
    'arm.pl.sup.nb1', {'3', '4', '5'};...
    'arm.pl.inf.nb1', {'arm.pl.sup.nb1'};...
    'arm.tr.sup.fiv', {'14:2:24'};...
    'arm.tr.inf.fiv', {'14:2:24'};...
    'arm.pl.inf.lock', {'''true''', '''false'''};...
    'arm.tr.inf.nb2',  {'1:4'};...
    'k_Minf.tr', {'1.4'};...
    'k_Minf.pl', {'1.4'};...
    'arm.tr.staffe.fi_sw', {'8'}
   %Compattabili
%     'L.X', {'34.69'};...
%     'L.Y', {'58.06'};...
%     'sollecitazioni.pl.q', {'50 - 7.45'};...
%     'sollecitazioni.pl.qSLU', {'1.3 * 8.8 + 1.5 * 50 - 30.23'};...
%     'num.campate' , {'8'};...
%     'num.pali_tr', {'13'};...
%     'arm.pl.sup.nb1', {'5', '5'};...
%     'arm.pl.inf.nb1', {'arm.pl.sup.nb1'};...
%     'arm.tr.sup.fiv', {'14:2:24'};...
%     'arm.tr.inf.fiv', {'14:2:24'};...
%     'arm.pl.inf.lock', {'''true''', '''false'''};...
%     'arm.tr.inf.nb2',  {'1:4'};...
%     'k_Minf.tr', {'2.3'};...
%     'k_Minf.pl', {'2.4'};...
%     'k_Msup.pl', {'1.8', '1.9', '2.0', '2.1' } 
%     'arm.tr.staffe.fi_sw', {'8', '10'}
    };
%%
risultati = table;
length_of_var2var = 1;
size_of_var2var = [];
head = regexprep(var2var(:,1), '\.', '_');
for i = 1:length(var2var)
    length_of_var2var = length_of_var2var * length(var2var{i, 2});
    size_of_var2var = [size_of_var2var, length(var2var{i, 2})];
    try
        evalin('base', ['risultati = [risultati, table({' var2var{i,1} '}, ''VariableNames'', {''' head{i} '''})];']);
    catch
        evalin('base', ['risultati = [risultati, table({' var2var{i,1} '}, ''VariableNames'', {''' head{i} '''})];']);
    end
end
risTab.ris = struct;
risTab = repmat(risTab, length_of_var2var, 1);
indice_min = [nan, nan];
L_min = [nan, nan];
% num_tr = nan;
% num_pali_tr = nan;
num_pali = nan;
costo_min = nan;
risultati = [risultati, table(costo_min, indice_min, L_min, num_pali)];
risultati = repmat(risultati, length_of_var2var, 1);
%%
main_waitbar = waitbar(0,'');
tempoTotale = 0;
try
    rmdir('Scenari', 's')
catch
    
end
    mkdir('Scenari')

for n_main = 1:length_of_var2var
    tic
    
    % aggiornamento waitbar
    waitbar(n_main/length_of_var2var, main_waitbar, sprintf('Scenario Corrente: %d/%d\nProgresso Globale: %.2f%%', n_main, length_of_var2var, n_main/length_of_var2var*100));
    set(main_waitbar, 'Name', sprintf('main: %.2f%%', n_main/length_of_var2var*100));
        
    % modifica delle variabili
    file2record = ['Scenari/scenario_' num2str(n_main) '.mat'];
    index_var2var = cell(1, length(var2var));
    [index_var2var{:}] = ind2sub(size_of_var2var, n_main);
    for i_var2var = 1:length(var2var)
        evalin('base', [var2var{i_var2var, 1} ' = ' var2var{i_var2var, 2}{index_var2var{i_var2var}} ';']);
        variable_ = evalin('base', var2var{i_var2var, 1});
        risultati{n_main, head{i_var2var}} = {variable_};
    end
   
    % calcolo
    core
    analisi
    save(file2record, var2save{:});
    
    risTab(n_main).ris = ris;
    risultati{n_main, 'indice_min'} = [i_min.pl, i_min.tr];
    risultati{n_main, 'L_min'} = [L.pl(i_min.pl), L.tr(i_min.tr)];
    risultati{n_main, 'num_pali'} = num.pali(i_min.pl, i_min.tr);    
    risultati{n_main, 'costo_min'} = costo_min
    save(file2record, 'M', 'V', 'soll', 'ris', '-append')
    
    % time passed
    elapsedTime = toc;
    tempoTotale = tempoTotale + elapsedTime;
    eta = tempoTotale/n_main * (length_of_var2var - n_main);
    timestr = sprintf('%-20s:\t%10.3f s\t(%s)\n', 'tempo iterazione', elapsedTime,  datestr(elapsedTime/(60*60*24),'HH:MM:SS'), 'tempo totale', tempoTotale, datestr(tempoTotale/(60*60*24),'HH:MM:SS'), 'eta', eta, datestr(eta/(60*60*24), 'HH:MM:SS'));
    disp(timestr)
end
risultati = [risultati, struct2table(risTab)];
%%

close(main_waitbar)
save('Scenari/sintesi','risultati')
%%
clearvars -except risultati
%%
[~, indexSorted] = sort(risultati.costo_min);
risSorted = risultati(indexSorted, :);
for i = 1:length(indexSorted)
    Scenari(i,1) = {sprintf('Scenari_%d', indexSorted(i))};
end
risSorted = [risSorted, table(Scenari)];
save('Scenari/sintesi', 'risSorted', '-append');
[costo_min, r] = min(risultati.costo_min);
load(['Scenari/Scenario_' num2str(r) '.mat']);

i_min.pl = ris.costo.ind_min(1);
i_min.tr = ris.costo.ind_min(2);
post