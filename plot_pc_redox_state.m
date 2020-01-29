function f = plot_pc_redox_state(species, ys, ts) % model specific variables)

pco = strcmp(species,'PCo');
pcr = strcmp(species,'PCr');

f = figure;
t = [];
PCo = [];
PCr = [];
for itime = 2:length(ys)
    PCo = [PCo, ys{itime}(pco,:)];
    PCr = [PCr, ys{itime}(pcr,:)];
    t = [t, ts{itime}];
end
plot(t, PCr./(PCo+PCr),'b')
title('PC pool redox state')
ylabel('reduced PC fraction')
xlabel('time')
    
    

