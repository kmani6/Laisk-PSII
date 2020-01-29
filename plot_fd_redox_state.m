function f = plot_fd_redox_state(species, ys, ts) % model specific variables)

fdo = strcmp(species,'FDo');
fdr = strcmp(species,'FDr');

f = figure;
hold on
t = [];
Y = [];
for itime = 2:length(ys)
    FDo = ys{itime}(fdo,:);
    FDr = ys{itime}(fdr,:);
    t= [t,ts{itime}];
    Y = [Y, FDr./(FDo+FDr)];
end
plot(t,Y);
title('Fd pool redox state')
ylabel('reduced Fd fraction')
xlabel('time')
    
    

