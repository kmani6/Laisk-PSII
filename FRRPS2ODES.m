function dydt = LaiskPS2ODES(t,y,k,k_init,rate_inds,S,Ynames,knames,parameters) 
% fprintf([num2str(t),'\n']);

n_flashes = parameters(1);
flash_duration = parameters(2);
flash_interval = parameters(3);
train_interval = parameters(4);
n_trains = parameters(5);
nrxn = length(rate_inds); 

% h = light_on(t);
r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = k(irxn)*prod(y(rate_inds{irxn}));
end
 
dydt = S*r; 

dFl = dLaiskFluorescence(Ynames,knames,k_init,y);
dydt(end+1) = dFl;

end


