
function [c,xopt] = TESTget_ms_parameters(x)

c = (x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2;
ceq = [];       
fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

lb = [0,0.2];
ub = [0.5,0.8];
A = [];
b = [];
Aeq = [];
beq = [];
x0 = [1/4,1/4];

nonlcon = @TESTget_ms_parameters;
                
ms = MultiStart('UseParallel', true);

model = @(x) fun(x);
 
rng default % For reproducibility
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon',...
                            'objective', model, ...
                            'x0', 'A', 'b', 'Aeq', 'beq', 'lb', 'ub', 'options', opts);
                        
                   
                                                      
[xopt] = run(ms,problem,3);


                

   
                