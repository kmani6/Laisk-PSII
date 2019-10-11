function GraphFig2(analyses_name)

            figure; 
            species_in_graph = {};
            colors = {'b','[0.75, 0.75, 0]','r'};
            
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        YrPrArBooIdx = find(contains(species, 'YrPrArBoo'));
        YoPrArBooIdx = find(contains(species, 'YoPrArBoo'));
        
        plot(t, Sol(YrPrArBooIdx,:),'o','color', colors{i},'MarkerSize',1.5);
        
        hold on
        plot(t, Sol(YoPrArBooIdx,:),'o','color', colors{i},'MarkerSize',1.5);
              
        species_in_graph{end+1} = [analyses_name{i},'YrPrArBoo'];
        species_in_graph{end+1} = [analyses_name{i}, 'YoPrArBoo'];
    end
    
        legend(species_in_graph);
        
end 