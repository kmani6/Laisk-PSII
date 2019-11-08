function GraphFig3(analyses_name)

        figure; 
            species_in_graph = {};
            colors = {'m','1, 0.5, 0]','r'};
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        [Fl1] = YoPrArBooLaiskFluorescence(species,knames,k,Sol);
        
         PSU1 = find(strcmp(knames,'PSU1'));
         a2 = find(strcmp(knames, 'a2'));
         Chl = find(strcmp(knames, 'Chl'));
         PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
                 
         semilogx(t(2:end), Fl(2:end),'o','color', colors{i},'MarkerSize',4);
         hold on
         semilogx(t(2:end), Fl1(2:end),'o','color',colors{i},'MarkerSize',4); 
       
%          species_in_graph{end+1} = [analyses_name{i},'Fl'];
%          species_in_graph{end+1} = [analyses_name{i}, 'Fl1'];
         
    end
    
         %legend({'jd 2000','jd 20','jd 200'}); 
         axis([.00001 .1 0 1.2])
%          ylim([0 1.2])
%          xlim([1e-5 .1])
         set(gca,'FontSize',22)
         set(gca,'linewidth',2)
         set(gca,'color','white')
                     
end 






