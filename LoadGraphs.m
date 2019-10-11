function [Graph] = LoadGraphDCMU(tspan,analyses_name)
        figure; 

    for i = 1:length(analyses_name) 
        load([analyses_name{i},'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
        [Fl] = LaiskFluorescence(Ynames,knames,k,Sol); 
        
          PSU1 = find(strcmp(knames,'PSU1'));
          a2 = find(strcmp(knames, 'a2'));
          Chl = find(strcmp(knames, 'Chl'));
          PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
          
          plot(t, Fl/PS2T);
          hold on

    end

    
    
end 
%         species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
%  
%         SumIndex1 = find(contains(species, 'YrPrAo'));
%         SumIndex2 = find(contains(species, 'YoPrAr'));
%         SumIndex3 = find(contains(species, 'YrPrAr'));
%         SumIndex4 = find(contains(species, 'YoPrAo'));
%         SumIndex5 = find(contains(species, 'YoPoAr'));
%         SumIndex6 = find(contains(species, 'YoPoAo'));
%         
%         semilogx(t, sum(Sol(SumIndex1,:)))
%         hold on
%         semilogx(t, sum(Sol(SumIndex2,:)))
%         hold on
%         semilogx(t, sum(Sol(SumIndex3,:)))
%         hold on
%         semilogx(t, sum(Sol(SumIndex4,:)))
%         hold on 
%         semilogx(t, sum(Sol(SumIndex5,:)))
%         hold on
%         semilogx(t, sum(Sol(SumIndex6,:)))
%         hold on
%         semilogx(t, Fl);
%         
%         legend(species_in_graph);

