function Print_Params(analysis_name)

file1 = [analysis_name,'/experimental_parameters'];
[parameters, parameter_names] = xlsread(file1);
n_flashes = parameters(1);
flash_duration = parameters(2);
flash_interval = parameters(3);
train_interval = parameters(4);
n_trains = parameters(5);

fprintf(n_flashes) = parameters(1); 

identifier = analysis_name(1:end);
fprintf(identifier,'\n',' n_flashes ',num2str(n_flashes), '\n',' flash_duration ',num2str(flash_duration),...
    '\n',' flash_interval ',num2str(train_interval),'\n',' n_trains ', num2str(n_trains))

end 



