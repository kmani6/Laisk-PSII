function [c,ceq] = nonlinconstr(x,kwfindep,kwbindep)
% Nonlinear equality constraints
ceq = x(kwfindep)/x(kwbindep) - 10^-14;
% Nonlinear inequality constraints
c = [];
end 

