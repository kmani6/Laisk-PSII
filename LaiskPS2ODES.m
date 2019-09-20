function dydt = LaiskPS2ODES(t,y,k,k_init,rate_inds,S,Ynames,knames) 
fprintf([num2str(t),'\n']);

nrxn = length(rate_inds); 

r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = k(irxn)*prod(y(rate_inds{irxn}));
end
 
dydt = S*r; 

dFl = dLaiskFluorescence(Ynames,knames,k_init,y);
dydt(end+1) = dFl;

end


