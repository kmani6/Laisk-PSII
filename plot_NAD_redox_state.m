function f = plot_NAD_redox_state(species, ys, ts) % model specific variables)

nadp = strcmp(species,'NADP');
nadph = strcmp(species,'NADPH');

f = figure;
hold on
NAD = [];
NADH = [];
t = [];
for itime = 2:length(ys)
    t = [t, ts{itime}];
    NAD = [NAD, ys{itime}(nadp,:)];
    NADH = [NADH, ys{itime}(nadph,:)];
end
plot(t, NADH./(NAD+NADH),'b')

title('NADPH pool redox state')
ylabel('reduced NAD fraction')
xlabel('time')
    
    

