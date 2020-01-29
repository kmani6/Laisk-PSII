function f = plot_pq_redox_state(species, ys, ts) % model specific variables)

pq = strcmp(species,'PQ');
pqh2 = strcmp(species,'PQH2');

f = figure;
t = [];
PQ = [];
PQH2 = [];

for itime = 2:length(ys)
    PQ = [PQ, ys{itime}(pq,:)];
    PQH2 = [PQH2, ys{itime}(pqh2,:)];
    t = [t, ts{itime}];
end
plot(t, PQH2./(PQ+PQH2))

title('PQ pool redox state')
ylabel('reduced PQ fraction')
xlabel('time')
    
    

