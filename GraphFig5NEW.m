function GraphFig5NEW(analyses_name)

        figure; 
        
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
        semilogx(t, Fl,'o','MarkerSize',3.5);
        hold on
         
    end
         
        legend({'Plastoquinone reduction red','Plastoquinone reduction blue'}); 
        axis([.00001 .01 0 1.2])
        set(gca,'color','white')
        set(gca,'FontSize',22)
        set(gca,'linewidth',2)


end 