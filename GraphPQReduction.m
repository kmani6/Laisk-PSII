function [Graph] = GraphFig2(tspan,analyses_name)
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