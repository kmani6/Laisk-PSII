function f = plot_H_species(species, ys, ts) % model specific variables)

hs = strcmp(species,'Hs');
hl = strcmp(species,'Hl');

f = figure;
hold on
t = [];
HS = [];
HL = [];
for itime = 2:length(ys)
    t = [t, ts{itime}];
    
    HS = [HS, ys{itime}(hs,:)];
    HL = [HL, ys{itime}(hl,:)];
    
end
plot(t, HL,'b')
plot(t, HS,'k')
legend({'H_{lumen}', 'H_{stroma}'})
title('Proton concentrations')
ylabel('H+ concentration')
xlabel('time')
    
    

