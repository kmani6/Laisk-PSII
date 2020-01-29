function f = plot_cyt_redox_state(species, ys, ts) % model specific variables)

Cyto = strcmp(species,'Cytfo');
Cytr = strcmp(species,'Cytfr');

f = figure;
% hold on
t = [];
C = [];
for itime = 2:length(ys)
    CYTo = ys{itime}(Cyto,:);
    CYTr = ys{itime}(Cytr,:);
    t = [t, ts{itime}];
    C = [C, CYTr./(CYTo+CYTr)];
   
end
plot(t,C);
title('Cyt pool redox state')
ylabel('reduced PC fraction')
xlabel('time')
    
    

