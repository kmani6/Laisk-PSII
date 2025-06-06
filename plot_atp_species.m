function f = plot_atp_species(species, ys, ts) % model specific variables)

atp = strcmp(species,'SynATP');
adp = strcmp(species,'SynADP');
F = strcmp(species,'SynF');
FD = strcmp(species,'SynFD');
FDP = strcmp(species,'SynFDP');
FT = strcmp(species,'SynFT');
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
legend({'ATP','ADP'})

f = figure;
hold on
for itime = 2:length(ys)
    F_all = ys{itime}(F,:);
    FD_all = ys{itime}(FD,:);
    FDP_all = ys{itime}(FDP,:);
    FT_all = ys{itime}(FT,:);
    
    plot(ts{itime}, F_all,'b')
    plot(ts{itime}, FD_all, 'k')
    plot(ts{itime}, FDP_all, 'm')
    plot(ts{itime}, FT_all, 'g')
end

title('ADP/ATP concentration')
ylabel('Concentration')
xlabel('time')
legend({'F', 'FD', 'FDP','FT'})
    

