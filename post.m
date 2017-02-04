%% Post processing
% plot dei grafici e scrittura dei risultati in file di testo
close all

try
    %% grafico delle frecce
    figure(1)
    subplot(1,2,1)
    surf(L.pl, L.tr, ris.tot.s_rel')
    title('Deformazioni massime')
    xlabel('L_{platea}')
    ylabel('L_{trave}')
    zlabel('(f_{pl}+f_{tr})/L_{tot}')
    hold on
    surf(L.pl, L.tr, f_lim*ones(length(L.pl), length(L.tr))')
    subplot(1,2,2)
    contourf(L.pl, L.tr, ris.tot.s_logico',[0, 1]);
    title('Soluzioni Valide')
    xlabel('L_{platea}')
    ylabel('L_{trave}')
    %% grafico dei costi
    figure(2)
    
    subplot(3,2,1)
    plot(L.pl, ris.costo.pl.tot)
    title('Costo Platea')
    xlabel('L_{platea}')
    ylabel('Costo Platea')
    
    subplot(3,2,3)
    [C.tr, h.tr] = contourf(L.pl, L.tr, ris.costo.tr.tot');
    clabel(C.tr, h.tr)
    title('Costo Travi')
    xlabel('L_{platea}')
    ylabel('L_{trave}')
    zlabel('Costo Travi')
    
    subplot(3,2,5)
    [C.pali, h.pali] = contourf(L.pl, L.tr, ris.costo.pali');
    clabel(C.pali, h.pali);
    title('Costo Pali')
    xlabel('L_{platea}')
    ylabel('L_{trave}')
    zlabel('Costo pali')
    
    subplot(3,2,[2, 4, 6])
    [C.tot, h.tot] = contourf(L.pl, L.tr, ris.costo.tot');
    clabel(C.tot, h.tot);
    title('Costo Totale')
    xlabel('L_{platea}')
    ylabel('L_{trave}')
    zlabel('Somma Costi')
catch
    disp('I risultati non possono essere mostrati graficamente')
end

%%
clc
out = fopen('results.txt','a+');
headfont(1:60) = {'#'};
fprintf(out, '%s', headfont{:});
fprintf(out, '\n\nFattore di amplificazione del momento per la platea:\n\tMinf x %.2f\n\tMsup x %.2f\n', k_Minf.pl, k_Msup.pl);
fprintf(out, '\n\nFattore di amplificazione del momento per la trave:\n\tMinf x %.2f\n\tMsup x %.2f\n', k_Minf.tr, k_Msup.tr);
fprintf(out, '\n\nLa soluzione ottimale si ha in corripondenza degli indici:\n');
fprintf(out, '\t% 4s\t% 8s\t% 8s\n', '', 'min', 'L[m]');
label = {'pl', 'tr'};
for i = 1:length(label)
    fprintf(out, '\ti_%2s\t% 8.0f\t% 8.2f \n', label{i}, i_min.(label{i}), L.(label{i})(i_min.(label{i})));
end
fprintf(out, '\nChe corrispondono a %d travi e %d pali per trave, per un totale di %d pali.\n', num.tr(i_min.pl), num.pali_tr(i_min.tr), num.pali(i_min.pl, i_min.tr));
fprintf(out, '\nLe deformazioni minime per questa combinazione sono pari a:\n');
fprintf(out, '\t% -4s\t% 12s\t% 12s\t% s<%.2e\n', '', 's_min[mm]', 's_rel[%]', 's_rel', f_lim );
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'pl', ris.pl.s_min(i_min.pl), ris.pl.s_rel(i_min.pl), ris.pl.s_logico(i_min.pl));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tr', ris.tr.s_min(i_min.pl, i_min.tr), ris.tr.s_rel(i_min.pl, i_min.tr), ris.tr.s_logico(i_min.pl, i_min.tr));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot1', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel(i_min.pl, i_min.tr), ris.tot.s_logico(i_min.pl, i_min.tr));
fprintf(out, '\t% -4s\t% 12.4f\t% 12.4e\t% 12d\n', 'tot2', ris.tot.s_min(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr), ris.tot.s_rel2(i_min.pl, i_min.tr)<= f_lim);
fprintf(out, '\nIl costo totale stimato è pari a %.2f €\n\n', ris.costo.tot(i_min.pl, i_min.tr));
fclose(out);
type('results.txt')