function[Graph] = MultipleGraphs(Var,analyses_name);

file = {};

% file = [analysis_name{i},'/LaiskReactions.xls'];

for i = 1:length(analyses_name)
    file{i} = dir(analyses_name{i});

    for ianalysis = 1:length(file)
    [~,reactions,~] = xlsread(file(ianalysis).name);
    [~,filename] = fileparts(files(ianalysis).name);
    
    save(filename, 'reactions')
    
    end 
    
end
 

