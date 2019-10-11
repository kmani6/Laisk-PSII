function GraphFig4(analyses_name)

        load([analyses_name,'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t');  
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
        species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
        
        PSU1 = find(strcmp(knames,'PSU1'));
        a2 = find(strcmp(knames, 'a2'));
        Chl = find(strcmp(knames, 'Chl'));
        PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
        
                 
        SumIndex1 = find(contains(species, 'YrPrAo'));
        SumIndex2 = find(contains(species, 'YoPrAr'));
        SumIndex3 = find(contains(species, 'YrPrAr'));
        SumIndex4 = find(contains(species, 'YoPrAo'));
        SumIndex5 = find(contains(species, 'YoPoAr'));
        SumIndex6 = find(contains(species, 'YoPoAo'));
      
        Color = 'magenta' ;
        S = scatter(t, sum(Sol(SumIndex1,:)/PS2T),[],Color);
        ylim([0 1.2])
        set(gca,'xscale','log')
        set(S, 'SizeData', 4)
        hold on
        
        Color = [1 .5 0];
        S = scatter(t, sum(Sol(SumIndex2,:)/PS2T),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 4)
        hold on
        
        Color = 'cyan' ;
        S = scatter(t, sum(Sol(SumIndex3,:)/PS2T),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 4)
        hold on
        
        Color = [.60 0 0] ;
        S = scatter(t, sum(Sol(SumIndex4,:)/PS2T),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 4)
        hold on
        
        Color = 'green' ;
        S = scatter(t, sum(Sol(SumIndex5,:)/PS2T),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 2)
        hold on
       
        Color = [.75 0 0] ;
        S = scatter(t, sum(Sol(SumIndex6,:)/PS2T),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 4)
        hold on
        
       legend(species_in_graph);
           
end 

