analysis_name = 'PFD15000_S_states';
file3 = [analysis_name,'/LaiskReactions.xls'];
[~,Rknames] = xlsread(file3);
new_reactions = {};
counter = 1;
reactions = Rknames(:,1);
nrxn = length(reactions);
ids = find(contains(reactions, 'Y'));
s_states = '01230';
species = {};
constants = {};
for irxn = 1:nrxn
    rstring = reactions{irxn}; %inds the current reaction in the rxns vector and sets as the rxn string
    if contains(reactions{irxn}, 'Y')
        sep_ind = strfind(rstring,'->'); %finds a string withing a string, and returns the idx
        lhs = strtrim(rstring(1:sep_ind-1)); %string trim b4 idx
        reactants = strtrim(strsplit(lhs,'+')); %splits the string into 2 components (sep. reactants, removes +)
        n_reactants = regexp(reactants,'(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
        rhs = strtrim(rstring(sep_ind+2:end));
        products = strtrim(strsplit(rhs,' + '));
        n_products = regexp(products, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
        if strcmp(Rknames{irxn,2}, 'jd')
            for i_state = 1:4
                tmp_reactants = {};
                for i_reactant = 1:length(n_reactants)
                    if contains(n_reactants{i_reactant}.species, 'Y')
                        tmp_reactant = [n_reactants{i_reactant}.stoich, ' ', 'S', s_states(i_state), n_reactants{i_reactant}.species];
                        tmp_reactants{end+1} = tmp_reactant;
                    else
                        tmp_reactant = [n_reactants{i_reactant}.stoich, ' ', n_reactants{i_reactant}.species];
                        tmp_reactants{end+1} = tmp_reactant;
                    end
                    if isempty(find(strcmp(species,tmp_reactant)))
                        species{end+1,1} = tmp_reactant;
                    end
                end
                tmp_reactants = strtrim(tmp_reactants);
                tmp_products = {};
                for i_product = 1:length(n_products)
                    if contains(n_products{i_product}.species, 'Y')
                        tmp_product = [n_products{i_product}.stoich, ' ', 'S', s_states(i_state+1), n_products{i_product}.species];
                        tmp_products{end+1} = tmp_product;
                        if isempty(find(strcmp(species,tmp_product)))
                            species{end+1,1} = tmp_product;
                        end
                    else
                        tmp_product = [n_products{i_product}.stoich, ' ', n_products{i_product}.species];
                        tmp_products{end+1} = tmp_product;
                        
                        if isempty(find(strcmp(species,tmp_product)))
                            species{end+1,1} = tmp_product;
                        end
                    end
                end
                tmp_products = strtrim(tmp_products);
                if i_state == 1
                    o2id = find(strcmp(tmp_products, 'O2'));
                    tmp_products(o2id) = [];
                    tmp_products{end+1} = 'Hl';
                elseif i_state == 2
                    o2id = find(strcmp(tmp_products, 'O2'));
                    tmp_products(o2id) = [];
                elseif i_state == 3
                    o2id = find(strcmp(tmp_products, 'O2'));
                    tmp_products(o2id) = [];
                    tmp_products{end+1} = 'Hl';
                elseif i_state == 4
                    tmp_products{end+1} = '2 Hl';
                end
                tmp_reactants = strjoin(tmp_reactants,' + ');
                tmp_products = strjoin(tmp_products,' + ');
                tmp_rxn = [tmp_reactants ' -> ' tmp_products];
                new_reactions{counter,1} = tmp_rxn;
                new_reactions{counter,2} = ['jd' s_states(i_state) s_states(i_state+1)];
                if isempty(find(strcmp(constants,new_reactions{counter,2})))
                    constants{end+1,1} = new_reactions{counter,2};
                end
                counter = counter+1;
            end
        else
            for i_state = 1:4
                tmp_reactants = {};
                for i_reactant = 1:length(n_reactants)
                    if contains(n_reactants{i_reactant}.species, 'Y')
                        tmp_reactant = [n_reactants{i_reactant}.stoich, ' ', 'S', s_states(i_state), n_reactants{i_reactant}.species];
                        tmp_reactants{end+1} = tmp_reactant;
                    else
                        tmp_reactant = [n_reactants{i_reactant}.stoich, ' ', n_reactants{i_reactant}.species];
                        tmp_reactants{end+1} = tmp_reactant;
                    end
                    if isempty(find(strcmp(species,tmp_reactant)))
                        species{end+1,1} = tmp_reactant;
                    end
                end
                tmp_reactants = strtrim(tmp_reactants);
                tmp_products = {};
                for i_product = 1:length(n_products)
                    if contains(n_products{i_product}.species, 'Y')
                        tmp_product = [n_products{i_product}.stoich, ' ', 'S', s_states(i_state), n_products{i_product}.species];
                        tmp_products{end+1} = tmp_product;
                    else
                        tmp_product = [n_products{i_product}.stoich, ' ', n_products{i_product}.species];
                        tmp_products{end+1} = tmp_product;
                        
                    end
                    if isempty(find(strcmp(species,tmp_product)))
                        species{end+1,1} = tmp_product;
                    end
                end
                tmp_products = strtrim(tmp_products);
                tmp_reactants = strjoin(tmp_reactants,' + ');
                tmp_products = strjoin(tmp_products,' + ');
                tmp_rxn = [tmp_reactants ' -> ' tmp_products];
                new_reactions{counter,1} = tmp_rxn;
                new_reactions{counter,2} = [Rknames{irxn,2}];
                if isempty(find(strcmp(constants,new_reactions{counter,2})))
                    constants{end+1,1} = new_reactions{counter,2};
                end
                counter = counter+1;
            end
            
        end
        
        
    else
        new_reactions{counter,1} = rstring;
        new_reactions{counter,2} = Rknames{irxn,2};
        if isempty(find(strcmp(constants,new_reactions{counter,2})))
            constants{end+1,1} = new_reactions{counter,2};
        end
        counter = counter+1;
        sep_ind = strfind(rstring,'->'); %finds a string withing a string, and returns the idx
        lhs = strtrim(rstring(1:sep_ind-1)); %string trim b4 idx
        reactants = strtrim(strsplit(lhs,'+')); %splits the string into 2 components (sep. reactants, removes +)
        n_reactants = regexp(reactants,'(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
        for i_reactant = 1:length(n_reactants)
            if isempty(find(strcmp(species,n_reactants{i_reactant}.species)))
                species{end+1,1} = n_reactants{i_reactant}.species;
            end
        end
        rhs = strtrim(rstring(sep_ind+2:end));
        products = strtrim(strsplit(rhs,' + '));
        n_products = regexp(products, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
        for i_product = 1:length(n_products)
            if isempty(find(strcmp(species,n_products{i_product}.species)))
                species{end+1,1} = n_products{i_product}.species;
            end
        end
    end
    
end

for irxn = 1:nrxn
    
    rstring = reactions{irxn}; %inds the current reaction in the rxns vector and sets as the rxn string
    sep_ind = strfind(rstring,'->'); %finds a string withing a string, and returns the idx
    lhs = strtrim(rstring(1:sep_ind-1)); %string trim b4 idx
    reactants = strtrim(strsplit(lhs,'+')); %splits the string into 2 components (sep. reactants, removes +)
    
    %WE SHOULD BE ABLE TO REMOVE STRTRIM IN 29
    
    n = regexp(reactants,'(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
    
    rate_inds{irxn} = []; %indices irxn in rate_inds and sets it as empty
    
    for ireactant = 1:length(n)
        
        if isempty(find(strcmp(species,n{ireactant}.species)))
            
            species{end+1} = n{ireactant}.species;
            
        end
        
        if ~isempty(n{ireactant}.stoich)
            s = regexp(n{ireactant}.stoich, '\(?([0-9.]+)\)?', 'tokens');
            stoich = str2double(s{1});
        else
            stoich = 1;
        end
        ireactant_ind = find(strcmp(species,n{ireactant}.species));
        Si = [Si;ireactant_ind];
        Sj = [Sj;irxn];
        Ss = [Ss; -stoich];
        rate_inds{irxn} = [rate_inds{irxn}; ireactant_ind];
        
    end
    
    
    rhs = strtrim(rstring(sep_ind+2:end));
    products = strtrim(strsplit(rhs,' + '));
    n = regexp(products, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
    
    for iproduct = 1:length(n)
        if ~isempty(n{iproduct})
            if isempty(find(strcmp(species,n{iproduct}.species)))
                species{end+1} = n{iproduct}.species;
            end
            if ~isempty(n{iproduct}.stoich)
                s = regexp(n{iproduct}.stoich, '\(?([0-9.]+)\)?', 'tokens');
                stoich = str2double(s{1});
            else
                stoich = 1;
            end
            
            iproduct_ind = find(strcmp(species,n{iproduct}.species));
            Si = [Si;iproduct_ind];
            Sj = [Sj;irxn];
            Ss = [Ss; stoich];
            
        end
    end
end