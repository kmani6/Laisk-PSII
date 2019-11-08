function GraphFig7(analyses_name)

        figure; 
            species_in_graph = {};
%             colors = {'r','b','g','cyan','brown','purple'};
            
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
        semilogx(t, Fl);
        hold on
         
    end
         
         legend({'Slow PSII','PFD600','PFD3000','PFD5000','PFD7500','PFD10000','PFD15000'}); 
         axis([.0001 1 0 1.2])
         set(gca,'color','white')
         set(gca,'FontSize',22)
         set(gca,'linewidth',2)
                     
end 




