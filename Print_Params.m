function Print_Params(analysis_name)

identifier = analysis_name(1:end);
fprintf('%s',identifier);

fprintf('analysis_name')

name = input('Enter your name: ');
first = name(1);
last = name(end);
fprintf('Your name is: %s %nWith first and last letters: %c and %c%n',name,first,last);



end 



