function GraphFig9(analyses_name)

        load([analyses_name,'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t');  
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
        species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
                 
        SumIndex1 = find(contains(species, 'YrPrAo'));
        SumIndex2 = find(contains(species, 'YoPrAr'));
        SumIndex3 = find(contains(species, 'YrPrAr'));
        SumIndex4 = find(contains(species, 'YoPrAo'));
        SumIndex5 = find(contains(species, 'YoPoAr'));
        SumIndex6 = find(contains(species, 'YoPoAo'));
      
        Color = 'magenta' ;
        S = scatter(t, sum(Sol(SumIndex1,:)),[],Color);
        ylim([0 1.2])
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = [1 .5 0];
        S = scatter(t, sum(Sol(SumIndex2,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = 'cyan' ;
        S = scatter(t, sum(Sol(SumIndex3,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = [.60 0 0] ;
        S = scatter(t, sum(Sol(SumIndex4,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = 'green' ;
        S = scatter(t, sum(Sol(SumIndex5,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
       
        Color = [.75 0 0] ;
        S = scatter(t, sum(Sol(SumIndex6,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 7)
        hold on
        
        semilogx(t, Fl,'o', 'color','red' ,'MarkerSize',4);
        hold on
        
         legend(species_in_graph);
         set(gca,'FontSize',20)
         set(gca,'linewidth',2)
         set(gca,'color','white')
         axis([.0001 .1 0 1])
            

end 