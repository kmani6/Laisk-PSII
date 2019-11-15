function Print_Params(analysis_name)

file1 = [analysis_name,'/experimental_parameters'];
[parameters, parameter_names] = xlsread(file1);
n_flashes = parameters(1);
flash_duration = parameters(2);
flash_interval = parameters(3);
train_interval = parameters(4);
n_trains = parameters(5);

identifier = analysis_name(1:end);
x = [newline,identifier,newline, '-----------------------------', newline,'n_flashes               ',num2str(n_flashes),newline,'flash_duration          ',num2str(flash_duration),...
    newline,'flash_interval          ',num2str(flash_interval),newline, 'train_interval          ',num2str(train_interval),newline,'n_trains                ',num2str(n_trains),newline]; 
disp(x)

end 



