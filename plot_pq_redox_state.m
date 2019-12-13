function f = plot_pq_redox_state(species, ys, ts) % model specific variables)

pq = strcmp(species,'PQ');
pqh2 = strcmp(species,'PQH2');

f = figure;
hold on
for itime = 2:length(ys)
    PQ = ys{itime}(pq,:);
    PQH2 = ys{itime}(pqh2,:);
    plot(ts{itime}, PQH2./(PQ+PQH2),'b')
end
title('PQ pool redox state')
ylabel('reduced PQ fraction')
xlabel('time')
    
    

