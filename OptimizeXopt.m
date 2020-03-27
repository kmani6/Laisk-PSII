function[fvalbestopt] = OptimizeXopt(analysis_name)
%Used as a "main" currently. Will make new directory to save everything in
%and will determine the best fval and xval. To change number of runs must 
%change the number in the while loop, as well as in Multistart in get_ms_paramaters. 

if ~exist('results', 'dir')
    mkdir('results')
end
if ~exist(['results/' analysis_name], 'dir')
    mkdir('results/', analysis_name);
end

datestring = datestr(now, 30); 

if ~exist(['results/' analysis_name '/' datestring], 'dir')
 	mkdir(['results/', analysis_name '/' datestring]);
end

    counter = 0;
    fvalbestopt = inf;
    while counter <= 2 
                  [xopt,fval] = get_ms_parameters(analysis_name, datestring);
                  
                    if fval < fvalbestopt 
                        fvalbestopt = fval;
                        xvalbestopt = xopt; 
                        counter = 0;
                    else 
                        counter = counter + 1; 
                    end 
          fprintf(counter)        
    end
    
    save(['results/', analysis_name '/' datestring '/optimal' datestr(now, 30)],...
    'fvalbestopt', 'xvalbestopt') 

    disp(fvalbestopt) 
    disp(xvalbestopt)
    
    
    
end


