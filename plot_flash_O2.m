function f = plot_S_states(species, ys, ts, Fs) % model specific variables)

o2 = contains(species,'O2');


flash_idcs = find(~isempty(Fs{i}));
O2f = zeros(length(flash_idcs));
f = figure;
hold on
for i = 1:length(flash_idcs)
    if ~isempty(Fs{i})
        O2f = 
    end
end
legend({'S0', 'S1', 'S2', 'S3'})
title('S-state distribution')
ylabel('S-state fraction')
xlabel('time, s')
    
    
    

