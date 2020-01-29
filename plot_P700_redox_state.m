function f = plot_P700_redox_state(species, ys, ts) % model specific variables)

P700o = strcmp(species,'P700o');
P700r = strcmp(species,'P700r');

f = figure;
% hold on
t = [];
C = [];
for itime = 2:length(ys)
    CYTo = ys{itime}(P700o,:);
    CYTr = ys{itime}(P700r,:);
    t = [t, ts{itime}];
    C = [C, CYTr./(CYTo+CYTr)];
end
plot(t,C);
title('P700 pool redox state')
ylabel('reduced P700 fraction')
xlabel('time')
    
    

