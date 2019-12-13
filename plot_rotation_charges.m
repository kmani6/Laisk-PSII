function f = plot_rotation_charges(species, ys, ts) % model specific variables)

r = strcmp(species,'rot');

f = figure;
hold on
for itime = 2:length(ys)
    R = ys{itime}(r,:);
    plot(ts{itime}, R,'b')
    
end
legend({'Rotation Charges'})
title('Rotation Charges')
ylabel('Rotation Charges')
xlabel('time')
    
    

