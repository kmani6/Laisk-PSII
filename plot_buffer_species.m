function f = plot_buffer_species(species, ys, ts) % model specific variables)
species_to_graph = {'HbufferS','bufferS', 'HbufferL','bufferL'};
% species_to_graph = {'Hs','Hl','Ha', 'Hb', 'ASYNHauuu',...
%                     'ASYNHbcuu','ASYNHacuu',...
%                     'ASYNHbccu','ASYNHaccu',...
%                     'ASYNHbccc'};
idcs = [];
for i = 1:length(species_to_graph)
    idcs(i) = find(strcmp(species,species_to_graph{i}));
end
colors = distinguishable_colors(length(species));


f = figure;
hold on
t = [];
Y = [];
for itime = 2:length(ys)
    t = [t, ts{itime}];
    Y = [Y, ys{itime}(idcs,:)];
end
plot(t, Y)
legend(species_to_graph)
title('buffer species')
ylabel('Buffer concentration')
xlabel('time')



