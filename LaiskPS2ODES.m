function dydt = LaiskPS2ODES(t,y,krxn,k,rate_inds,S,Ynames,knames)

% if (y(PQH2)~=0)
%     krxn(rqr1) = krxn(rqr1)*y(PQH2);
%     
% end
% 
% if (y(PQ)~=0)
%     krxn(oqr1) = krxn(oqr1)*y(PQ);
%     
% end


% fprintf([num2str(t),'\n']);

nrxn = length(rate_inds);

r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end

dydt = S*r;

dFl = dLaiskFluorescence(Ynames,knames,k,y);
dydt(end+1) = dFl;

end


