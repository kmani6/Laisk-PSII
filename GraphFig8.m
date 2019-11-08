function GraphFig8(analyses_name)

        figure; 
    
        load(['K-peak','/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        SumIndex1 = find(contains(species, 'YoPrAr'));
        SumIndex2 = find(contains(species, 'YoPrAo'));
        SumIndex3 = find(contains(species, 'YoPoAo'));
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
        Color = [1 .5 0];
        S = scatter(t, sum(Sol(SumIndex1,:)),[],Color);
        ylim([0 1.2])
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = 'magenta';
        S = scatter(t, sum(Sol(SumIndex2,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
        
        Color = 'cyan' ;
        S = scatter(t, sum(Sol(SumIndex3,:)),[],Color);
        set(gca,'xscale','log')
        set(S, 'SizeData', 8)
        hold on
       
        semilogx(t, Fl,'o', 'color','red' ,'MarkerSize',4);
        
        legend({'YoPrAr','YoPrAo','YoPoAo','Fl'}); 
        axis([.0001 .2 0 1.2])
        set(gca,'color','white')
        set(gca,'FontSize',22)
        set(gca,'linewidth',2)
        hold on
        
        species_in_graph = {{'PQH2','PQ'}};
        lgd = {};
        
        for i = 1:length(species_in_graph)
    
            current_species = species_in_graph{i};
            idcs = [];
                for j = 1:length(current_species)
                    idcs(end+1) = find(strcmp(species, current_species{j}));        
                end
                
        redox_state = Sol(idcs(1),:)./(Sol(idcs(1),:) + Sol(idcs(2),:));
        plot(t,redox_state,'o','MarkerSize',3.5);
        lgd{end+1} = (strcat(current_species{1}, "/(", current_species{2}, " + ", current_species{1}, ")"))  ;
        
        end
        
        plot(t,Fl,'o','MarkerSize',3.5);
        set(gca,'FontSize',20)
        set(gca,'linewidth',2)
        ylim([0 1.2])

end 
