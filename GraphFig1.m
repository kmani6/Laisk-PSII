function GraphFig1(analyses_name)

        figure; 
            species_in_graph = {};
            colors = {'b','0.75, 0.75, 0]','r'};
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        [Fl1] = YoPrArBooLaiskFluorescence(species,knames,k,Sol);
        
         PSU1 = find(strcmp(knames,'PSU1'));
         a2 = find(strcmp(knames, 'a2'));
         Chl = find(strcmp(knames, 'Chl'));
         PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
                 
         plot(t, Fl,'o', 'color', colors{i},'MarkerSize',2);
         hold on
         plot(t, Fl1,'o','color',colors{i},'MarkerSize',2); 
         

         
         species_in_graph{end+1} = [analyses_name{i},'Fl'];
         species_in_graph{end+1} = [analyses_name{i}, 'Fl1'];
         
    end
    
         legend(species_in_graph); 
         ylim([0 1.2])
         xlim([0 .0002])
%          set(gca,'color','none')
                     
end 


