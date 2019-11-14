[Sol] = FRRFluorescence2(analysis_name)
    
    load([analyses_name,'/result.mat'],'Sol','k','knames','y0','Ynames','species','kconst','t');  
    
    PFD = find(strcmp(knames, 'PFD')); 
    Labs = find(strcmp(knames,'Labs'));
    a2 = find(strcmp(knames, 'a2'));
    PSU1 = find(strcmp(knames,'PSU1'));
    Chl = find(strcmp(knames, 'Chl')); 
  
    PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
    n1 = k(PFD)*k(Labs)*(1-k(a2))/PS1T;
    n2 = k(PFD)*k(Labs)*k(a2)/PS2T; 
    
    mult1Val = n2*k(kp)/(1+k(kp)+k(kn)+k(kr));
    mult2Val = n2*k(kp)/(1+k(kp)+k(kn));
    
    k(mult1) = mult1Val;
    k(mult2) = mult2Val;
    k(n1id) = n1; 
    
    ys = {};
    Fs = {};
    ts = {};
    PQ = find(strcmp(species,'PQ'));
    PQH2 = find(strcmp(species, 'PQH2'));
    
    % dark adapt the system
    k(mult1) = 0;
    k(mult2) = 0;    
    k(n1id) = 0;

        %         t = logspace(-8, log10(flash_duration), 200);
        %         t(1) = 0;
dark_adaptation_time = 180;
t = linspace(0, dark_adaptation_time, dark_adaptation_time*1000);
tic;
Sol =  ode2(@(t,y) FRRPS2ODES(t,y,k(kconst),k,kconst,rate_inds,S,species,knames,PQ,PQH2,oqrinds,rqrinds,species,Rknames(:,1)),t,yinitial);
toc
Sol = Sol';
ys{end+1} = Sol; %calculate the species evolution during the light
ts{end+1} = -dark_adaptation_time + t; %save the times used
Fs{end+1} = LaiskFluorescence(species,knames,k,Sol); %save the fluoprescence values
yinitial = Sol(:,end); %initialize the y vector for the next iteration 

% for i =1:size(Sol,1)
% figure;
% plot(t,Sol(i,:));
% legend(species(i));
% end

kmod = k(kconst);
nrxn = length(Rknames);

r = [];
ttmp = cell2mat(ts(end));
ytmp = cell2mat(ys(end));
for irxn = 1:nrxn
    for j =1:length(ttmp)
        r(irxn,j) = kmod(irxn)*prod(ytmp(rate_inds{irxn},j));
    end
    figure; 
    subplot(length(rate_inds{irxn})+1,1,1)
    plot(ttmp, r(irxn,:))
    legend(Rknames(irxn,1))
    for k = 1:length(rate_inds{irxn})
        subplot(length(rate_inds{irxn})+1,1,k+1)
        plot(ttmp, ytmp(rate_inds{irxn}(k),:))
        legend(species(rate_inds{irxn}(k)))
    end
end


end 