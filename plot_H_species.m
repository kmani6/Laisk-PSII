function f = plot_H_species(species, ys, ts) % model specific variables)

hs = strcmp(species,'Hs');
hl = strcmp(species,'Hl');

f = figure;
hold on
for itime = 2:length(ys)
    HS = ys{itime}(hs,:);
    HL = ys{itime}(hl,:);
    plot(ts{itime}, HL,'b')
    plot(ts{itime}, HS,'k')
    
end
legend({'H_{lumen}', 'H_{stroma}'})
title('Proton concentrations')
ylabel('H+ concentration')
xlabel('time')
    
    

