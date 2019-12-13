function f = plot_NAD_redox_state(species, ys, ts) % model specific variables)

nadp = strcmp(species,'NADP');
nadph = strcmp(species,'NADPH');

f = figure;
hold on
for itime = 2:length(ys)
    NAD = ys{itime}(nadp,:);
    NADH = ys{itime}(nadph,:);
    plot(ts{itime}, NADH./(NAD+NADH),'b')
end
title('NADH pool redox state')
ylabel('reduced NAD fraction')
xlabel('time')
    
    

