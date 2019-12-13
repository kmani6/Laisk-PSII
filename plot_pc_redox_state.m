function f = plot_pc_redox_state(species, ys, ts) % model specific variables)

pco = strcmp(species,'PCo');
pcr = strcmp(species,'PCr');

f = figure;
hold on
for itime = 2:length(ys)
    PCo = ys{itime}(pco,:);
    PCr = ys{itime}(pcr,:);
    plot(ts{itime}, PCr./(PCo+PCr),'b')
end
title('PC pool redox state')
ylabel('reduced PC fraction')
xlabel('time')
    
    

