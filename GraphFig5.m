function GraphFig5(analyses_name)

        figure; 
            species_in_graph = {};
            colors = {'b','r'};
            
    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        
        [Fl] = LaiskFluorescence(species,knames,k,Sol); 
        
         PSU1 = find(strcmp(knames,'PSU1'));
         a2 = find(strcmp(knames, 'a2'));
         Chl = find(strcmp(knames, 'Chl'));
         PQT = find(strcmp(knames,'PQT'));
         PQH2 = find(strcmp(species, 'PQH2'));
         PS2T = (1-k(a2))*(k(Chl)/k(PSU1)); 
        
%         PQredstate = Sol(PQH2)/Sol(PQ);
         
         if (k(PQT)~=0)
             
            PQredstate = Sol(PQH2)/k(PQT);
            plot(PQredstate, Fl/PS2T,'o', 'color', colors{i},'MarkerSize',3.5);
            hold on

         end  
               
%         plot(PQredstate, Fl/PS2T,'o', 'color', colors{i},'MarkerSize',3.5);
%         hold on
%         
    end
         
         legend({'rqr100,oqd100','rqr2000,oqd2000'}); 
         ylim([0 1.2])
         xlim([0 1])
         set(gca,'color','white')
                     
end 






