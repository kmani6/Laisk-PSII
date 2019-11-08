function GraphFig6(analyses_name)

load(['full_induction_test','/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t')  
[Fl] = LaiskFluorescence(species,knames,k,Sol); 
 
figure;
hold on
species_in_graph = {{'PQH2','PQ'}, {'Cytfr','Cytfo'} {'PCr','PCo'}, {'P700r','P700o'}, {'FDr','FDo'}};
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
    
plot(t,Fl,'o','MarkerSize',3.5);
set(gca,'FontSize',20)
set(gca,'linewidth',2)
ylim([0 1.2])
lgd{end+1} = "Fl";
legend(lgd);
 
end 