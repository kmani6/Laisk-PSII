function r = LaiskRates(t,y,k,rate_inds,S) 
 
nrxn = length(rate_inds); 

r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = k(irxn)*prod(y(rate_inds{irxn}));
    if irxn == 25
        foo = 1;
    end
end
end


