function dydt = PS2ODES(y,krxn,k,rate_inds,S,Ynames,knames)

nrxn = length(rate_inds);

r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end

dydt = S*r;

dFl = dLaiskFluorescence(Ynames,knames,k,y);
dydt(end+1) = dFl;

end


