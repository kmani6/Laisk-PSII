function f = plot_atp_species(species, ys, ts) % model specific variables)

atp = strcmp(species,'SynATP');
adp = strcmp(species,'SynADP');

f = figure;
hold on
for itime = 2:length(ys)
    ATP = ys{itime}(atp,:);
    ADP = ys{itime}(adp,:);
    plot(ts{itime}, ATP./(ATP+ADP),'b')
end
title('ADP phosphorylation state')
ylabel('ATP fraction')
xlabel('time')

f = figure;
hold on
for itime = 2:length(ys)
    ATP = ys{itime}(atp,:);
    ADP = ys{itime}(adp,:);
    plot(ts{itime}, ATP,'b')
    plot(ts{itime}, ADP, 'k')
end
title('ADP/ATP concentration')
ylabel('Concentration')
xlabel('time')
legend({'ATP', 'ADP'})
    

