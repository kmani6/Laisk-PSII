function dydt = LaiskPS2ODES(t,y,k,k_init, kinds,rate_inds,S,Ynames,knames,PQ,PQH2,oqrinds,rqrinds,species, reactions) 
% fprintf([num2str(t),'\n']);

nrxn = length(rate_inds); 

k(oqrinds) = k(oqrinds)/y(PQ);
k(rqrinds) = k(rqrinds)*y(PQH2);

% h = light_on(t);
r = zeros(nrxn,1);

for irxn = 1:nrxn
    r(irxn,1) = k(irxn)*prod(y(rate_inds{irxn}));
end
 
dydt = S*r; 
% reactions_to_check = [21,27];
% for i = reactions_to_check
% %     disp(i)
%     disp(reactions{i})
%     disp(strcat('   rate = ', num2str(r(i))))
%     for j = 1:length(rate_inds{i})
%         disp(strcat(species(rate_inds{i}(j)),'=',num2str(y(rate_inds{i}(j)))))
%     end
%     disp(strcat(knames(kinds(i)), ' = ', num2str(k(i))))
% end

end


