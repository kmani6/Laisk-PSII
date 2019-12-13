function dydt = PS2ODES(t,y,krxn,k,rate_inds,S,Rknames,species)

nrxn = length(rate_inds);

r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end

dydt = S*r;

% dFl = 0;% dLaiskFluorescence(Ynames,knames,k,y);
% dydt(end+1) = dFl;

end


