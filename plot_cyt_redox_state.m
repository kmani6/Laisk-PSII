function f = plot_cyt_redox_state(species, ys, ts) % model specific variables)

Cyto = strcmp(species,'Cytfo');
Cytr = strcmp(species,'Cytfr');

f = figure;
hold on
for itime = 2:length(ys)
    CYTo = ys{itime}(Cyto,:);
    CYTr = ys{itime}(Cytr,:);
    plot(ts{itime}, CYTr./(CYTo+CYTr),'b')
end
title('Cyt pool redox state')
ylabel('reduced PC fraction')
xlabel('time')
    
    

