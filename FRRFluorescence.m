function FRRFluorescence(analysis_name)

file1 = [analysis_name,'/experimental_parameters'];
[parameters, parameter_names] = xlsread(file1);
n_flashes = parameters(1);
flash_duration = parameters(2);
flash_interval = parameters(3);
train_interval = parameters(4);
n_trains = parameters(5);

maxtime = n_trains*(n_flashes*(flash_duration+ flash_interval) + train_interval);
 
file1 = [analysis_name,'/LaiskConstants.xls'];
[k,knames] = xlsread(file1);
 
file2 = [analysis_name,'/LaiskY.xls'];
[y0,Ynames] = xlsread(file2);
 
file3 = [analysis_name,'/LaiskReactions.xls'];
[~,Rknames] = xlsread(file3);
 
%k = rand(37,1);
%y0 = rand(38,1); 
 
% a2 = find(strcmp(knames,'a2'));
% Chl = find(strcmp(knames,'Chl'));
% PSU1 = find(strcmp(knames,'PSU1'));
% P700o = find(strcmp(Ynames, 'P700o'));
% FXo = find(strcmp(Ynames, 'FXo'));
% 
% PS1T = (1-k(a2))*(k(Chl)/k(PSU1));
% y(P700o) = PS1T;
% y(FXo) = PS1T; 
 
PFD = find(strcmp(knames, 'PFD')); 
a2 = find(strcmp(knames, 'a2'));
b1d = find(strcmp(knames, 'b1d'));
b2d = find(strcmp(knames, 'b2d'));
Chl = find(strcmp(knames, 'Chl')); 
CytfT = find(strcmp(knames, 'CytfT')); 
FDT = find(strcmp(knames, 'FDT')); 
jd = find(strcmp(knames, 'jd')); 
kb6f = find(strcmp(knames, 'kb6f')); 
kcytf = find(strcmp(knames, 'kcytf')); 
kf = find(strcmp(knames, 'kf')); 
kfd = find(strcmp(knames, 'kfd')); 
kfx = find(strcmp(knames, 'kfx')); 
kn = find(strcmp(knames, 'kn')); 
kp = find(strcmp(knames, 'kp'));
kpc = find(strcmp(knames, 'kpc'));
kr = find(strcmp(knames, 'kr')); 
kE1 = find(strcmp(knames, 'kE1'));
kEb6f = find(strcmp(knames, 'kEb6f'));
kEcytf = find(strcmp(knames,'kEcytf'));
kEfx = find(strcmp(knames,'kEfx'));
kEpc = find(strcmp(knames,'kEpc'));
Labs = find(strcmp(knames,'Labs'));
oqd = find(strcmp(knames,'oqd'));
oqr = find(strcmp(knames,'oqr'));
PCT = find(strcmp(knames,'PCT'));
PQT = find(strcmp(knames,'PQT'));
PSU1 = find(strcmp(knames,'PSU1'));
PSU2 = find(strcmp(knames,'PSU2'));
rqd = find(strcmp(knames,'rqd'));
rqr = find(strcmp(knames,'rqr'));
b1r = find(strcmp(knames,'b1r'));
kq = find(strcmp(knames,'kq')); 
P700T = find(strcmp(knames, 'P700T'));
FXT = find(strcmp(knames, 'FXT')); 
 
YoPoAo = find(strcmp(Ynames,'YoPoAo')); 
YoPoAoBoo = find(strcmp(Ynames,'YoPoAoBoo')); 
YoPrAo = find(strcmp(Ynames,'YoPrAo')); 
YoPoAr = find(strcmp(Ynames,'YoPoAr')); 
YoPrAoBoo = find(strcmp(Ynames,'YoPrAoBoo')); 
YoPoArBoo = find(strcmp(Ynames,'YoPoArBoo')); 
YoPoAoBro = find(strcmp(Ynames,'YoPoAoBro')); 
YrPrAo = find(strcmp(Ynames,'YrPrAo')); 
YoPrAr = find(strcmp(Ynames,'YoPrAr')); 
YrPrAoBoo = find(strcmp(Ynames,'YrPrAoBoo')); 
YoPrArBoo = find(strcmp(Ynames,'YoPrArBoo')); 
YoPrAoBro = find(strcmp(Ynames,'YoPrAoBro')); 
YoPoArBro = find(strcmp(Ynames,'YoPoArBro')); 
YoPoAoBrr = find(strcmp(Ynames,'YoPoAoBrr')); 
YrPrAr = find(strcmp(Ynames,'YrPrAr')); 
YrPrArBoo = find(strcmp(Ynames,'YrPrArBoo')); 
YoPrArBro = find(strcmp(Ynames,'YoPrArBro')); 
YoPrAoBrr = find(strcmp(Ynames,'YoPrAoBrr')); 
YrPrAoBro = find(strcmp(Ynames,'YrPrAoBro')); 
YoPoArBrr = find(strcmp(Ynames,'YoPoArBrr')); 
YrPrArBro = find(strcmp(Ynames,'YrPrArBro')); 
YrPrAoBrr = find(strcmp(Ynames,'YrPrAoBrr')); 
YoPrArBrr = find(strcmp(Ynames,'YoPrArBrr')); 
YrPrArBrr = find(strcmp(Ynames,'YrPrArBrr')); 
PQH2 = find(strcmp(Ynames,'PQH2')); 
PQ = find(strcmp(Ynames,'PQ')); 
Cytfr = find(strcmp(Ynames,'Cytfr')); 
PCr = find(strcmp(Ynames,'PCr')); 
P700r = find(strcmp(Ynames,'P700r')); 
FXr = find(strcmp(Ynames,'FXr')); 
FDr = find(strcmp(Ynames,'FDr')); 
PCo = find(strcmp(Ynames, 'PCo'));
P700o = find(strcmp(Ynames, 'P700o'));
FXo = find(strcmp(Ynames, 'FXo'));
FDo = find(strcmp(Ynames, 'FDo'));
Cytfo = find(strcmp(Ynames, 'Cytfo')); 
O2 = find(strcmp(Ynames, 'O2')); 
Fl = find(strcmp(Ynames, 'Fl')); 
 
PS1T = k(a2)*k(Chl)/k(PSU2); 
PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
PS2TR = (1-k(a2))*(k(Chl)/k(PSU1));

% PS1T = 1; 
% PS2T = 1;
n1 = k(PFD)*k(Labs)*(1-k(a2))/PS1T;
n2 = k(PFD)*k(Labs)*k(a2)/PS2T; 
y0(P700r) = PS1T/PS2T;
y0(FXo) = PS1T/PS2T;


k(oqr) = k(oqr);
k(rqr) = k(rqr);
k(kb6f) = k(kb6f);
k(kcytf) = k(kcytf);
k(kpc) = k(kpc);
k(kfx) = k(kfx);
k(P700T) = PS1T/PS2T; 
k(FXT) = PS1T/PS2T; 
k(kfd) = k(kfd);

 
%Vb6f = k(kb6f)*(y0(PQH2)*y0(Cytfo)-y0(Cytfr)*(y0(PQ)/k(kEb6f)));
 
mult1 = find(strcmp(knames,'n2*kp/(1+kp+kn+kr)'));
mult2 = find(strcmp(knames,'n2*kp/(1+kp+kn)')); 
Div1 = find(strcmp(knames,'kpc/kEpc'));
Div2 = find(strcmp(knames,'kcytf/kEcytf'));
Div3 = find(strcmp(knames,'kfx/kEfx'));
Div4 = find(strcmp(knames,'kb6f/kEb6f'));
n1idx = find(strcmp(knames,'n1'));

mult1Val = n2*k(kp)/(1+k(kp)+k(kn)+k(kr));
mult2Val = n2*k(kp)/(1+k(kp)+k(kn));
Div1Val = k(kpc)/k(kEpc);
Div2Val = k(kcytf)/k(kEcytf);
Div3Val = k(kfx)/k(kEfx);
Div4Val = k(kb6f)/k(kEb6f);
 
k(mult1) = mult1Val;
k(mult2) = mult2Val;
k(Div1) = Div1Val;
k(Div2) = Div2Val;
k(Div3) = Div3Val; 
k(Div4) = Div4Val;
k(n1idx) = n1;


 
[species,S,rate_inds] = Laisk_read_excel_model(analysis_name);
yinitial = zeros(length(y0),1);
for i = 1:length(Ynames)
    index = find(strcmp(species,Ynames(i)));
    yinitial(index) = y0(i);
end

[kconst] = LaiskKconstants(analysis_name);
 
%THE FOLLOWING CODE NEEDS TO GO IN ITS OWN FUNCTION IN ORDER TO KEEP THINGS
%CONCISE-------------------------------------------------------------------
ys = {};
Fs = {};
ts = {};
PQ = find(strcmp(species,'PQ'));
PQH2 = find(strcmp(species, 'PQH2'));

% dark adapt the system
k(mult1) = 0;
k(mult2) = 0;    
k(n1idx) = 0;
%t = logspace(-8, log10(flash_duration), 200);
%t(1) = 0;
dark_adaptation_time = 30; %3-5 minutes typically
t = linspace(0, dark_adaptation_time, dark_adaptation_time*1e5);
tic;
Sol =  ode2(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
toc
Sol = Sol';
ys{end+1} = Sol; %calculate the species evolution during the light
ts{end+1} = -dark_adaptation_time + t; %save the times used
Fs{end+1} = []; %save an empty fluorescence vector
yinitial = Sol(:,end); %initialize the y vector for the next iteration 

% for i =1:size(Sol,1)
% figure;
% plot(t,Sol(i,:));
% legend(species(i));
% end
% 
% kmod = k(kconst);
% nrxn = length(Rknames);
% 
% r = [];
% ttmp = cell2mat(ts(end));
% ytmp = cell2mat(ys(end));
% for irxn = 1:nrxn
%     for j =1:length(ttmp)
%         r(irxn,j) = kmod(irxn)*prod(ytmp(rate_inds{irxn},j));
%     end
%     figure; 
%     subplot(length(rate_inds{irxn})+1,1,1)
%     plot(ttmp, r(irxn,:))
%     legend(Rknames(irxn,1))
%     for k = 1:length(rate_inds{irxn})
%         subplot(length(rate_inds{irxn})+1,1,k+1)
%         plot(ttmp, ytmp(rate_inds{irxn}(k),:))
%         legend(species(rate_inds{irxn}(k)))
%     end
% end
% 

for train = 1:n_trains
    fprintf('train %i \n', train)
    for flash = 1:n_flashes
        fprintf('flash %i \n', flash)

        k(mult1) = mult1Val;
        k(mult2) = mult2Val;    
        k(n1idx) = n1;
%         t = logspace(-8, log10(flash_duration), 200);
%         t(1) = 0;
        t = linspace(0, flash_duration, flash_duration*1e6);
        tic;
        Sol = ode2(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
        toc
        Sol = Sol';
        ys{end+1} = Sol; %calculate the species evolution during the light
        ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + (flash-1)*(flash_duration+flash_interval)+t; %save the times used
        Fs{end+1} = LaiskFluorescence(species,knames,k,Sol); %save the fluoprescence values
        yinitial = Sol(:,end); %initialize the y vector for the next iteration 
        
        k(mult1) = 0;
        k(mult2) = 0;      
        k(n1idx) = 0;
%         t = logspace(-5, log10(flash_interval), 200); %assign the time interval appropriate for the dark interval)
%         t(1) = 0;
        t = linspace(0, flash_interval, flash_interval*1e6);
        tic;
        Sol = ode2(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
        toc
        tic;
        Sol = Sol';
        toc
        ys{end+1} = Sol; %calculate the species evolution during the dark between flashes
        ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + (flash-1)*(flash_duration+flash_interval) + flash_duration +t; %save the times used
        Fs{end+1} = []; % save an empty vector for fluorescence (to allign the values of times and fluorescence)
        yinitial = Sol(:,end); %initialize the y vector for the next iteration 
%         subplot(1,3,1); plot(ts{end-1},Fs{end-1})
%         foo = ys{end-1}; subplot(1,3,2); plot(ts{end-1},foo(1:end-1,:))
%         foo = ys{end}; subplot(1,3,3); plot(ts{end},foo(1:end-1,:))
        toc
        pause(eps)
    end
    k(mult1) = 0;
    k(mult2) = 0;       
    k(n1idx) = 0;
%     t = logspace(-5, log10(0.025), 100); %assign the time interval appropriate for the dark interval)
%     t(1) = 0;
    t = linspace(0, train_interval, train_interval*1e5);
    tic;
    Sol = ode2(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
    toc
    Sol = Sol';
    ys{end+1} = Sol; %calculate the species evolution during the dark between flashes
    ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+t; %save the times used
    Fs{end+1} = []; % save an empty vector for fluorescence (to allign the values of times and fluorescence)
%     SumIndex1 = find(contains(species, 'YrPrAo'));
%     SumIndex2 = find(contains(species, 'YoPrAr'));
%     SumIndex3 = find(contains(species, 'YrPrAr'));
%     SumIndex4 = find(contains(species, 'YoPrAo'));
%     SumIndex5 = find(contains(species, 'YoPoAr'));
%     SumIndex6 = find(contains(species, 'YoPoAo'));
%     figure;
%     species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo'};
%     plot(ts{end}, sum(Sol(SumIndex1,:)),'r')
%     hold on
%     plot(ts{end}, sum(Sol(SumIndex2,:)),'g')
%     plot(ts{end}, sum(Sol(SumIndex3,:)),'b')
%     plot(ts{end}, sum(Sol(SumIndex4,:)),'c')
%     plot(ts{end}, sum(Sol(SumIndex5,:)),'m')
%     plot(ts{end}, sum(Sol(SumIndex6,:)),'k')
%     hold off
%     legend(species_in_graph); 
    yinitial = Sol(:,end); %initialize the y vector for the next iteration    
end

save([analysis_name '/FRR_results.mat'], 'ts', 'Fs', 'ys', 'knames', 'kconst', 'k', 'species', 'Rknames','-v7.3')  
% 
% figure;
% hold on
% for i = 1:length(ts)
%     if length(ts{i}) == length(Fs{i})
%         plot(ts{i}, Fs{i}, 'r')
%     end
% end
% legend('Fluorescence');


% t = logspace(-5, -2, 1000);
% t(1) = 0;
% Sol = ode2(@(t,y) FRRPS2ODES(t,y,k(kconst),k,rate_inds,S,species,knames,parameters),t,yinitial);
% dydt = [];
% for i = 1:length(Sol.x)
%     dydt(:,i) = LaiskPS2ODES(Sol.x(i),Sol.y(:,i),k(kconst),rate_inds,S);
%     r(:,i) = LaiskRates(Sol.x(i),Sol.y(:,i),k(kconst),rate_inds,S);
% end


% graph_colors = 'bgrcmyk';
% figure;
% species_in_graph = {'PQ', 'PQH2'};
% hold on
% idcs = [];
% for i = 1:length(species_in_graph)
%     idcs(end+1) = find(strcmp(species, species_in_graph(i)));
% end
% for i = 1:length(ts)
%     y = ys{i};
%     for j = 1:length(species_in_graph)
%         plot(ts{i},y(idcs(j),:),graph_colors(j))
%     end
% end
% legend(species_in_graph)
% 
% for s = 1:length(species)
%     if rem(s-1,4) == 0
%         figure;
%     end
%     if rem(s,4) == 0
%         pos = 4;
%     else
%         pos = rem(s,4);
%     end
%     subplot(2,2,pos)
%     hold on
%     for i = 1:length(ts)
%         y = ys{i};
%         plot(ts{i},y(s,:),'k')
%     end
%     legend(species(s));
% end
nrxn = length(Rknames);
% kmod = k(kconst);
% r = [];
% ttmp = cell2mat(ts(end));
% ytmp = cell2mat(ys(end));
% % for irxn = 1:nrxn
%     for j =1:length(ttmp)
%         r(irxn,j) = kmod(irxn)*prod(ytmp(rate_inds{irxn},j));
%     end
%     figure; 
%     subplot(length(rate_inds{irxn})+1,1,1)
%     plot(ttmp, r(irxn,:))
%     legend(Rknames(irxn,1))
%     for k = 1:length(rate_inds{irxn})
%         subplot(length(rate_inds{irxn})+1,1,k+1)
%         plot(ttmp, ytmp(rate_inds{irxn}(k),:))
%         legend(species(rate_inds{irxn}(k)))
%     end
% end

 SumIndex1 = find(contains(species, 'YrPrAo'));
 SumIndex2 = find(contains(species, 'YoPrAr'));
 SumIndex3 = find(contains(species, 'YrPrAr'));
 SumIndex4 = find(contains(species, 'YoPrAo'));
 SumIndex5 = find(contains(species, 'YoPoAr'));
 SumIndex6 = find(contains(species, 'YoPoAo'));
%  
% figure; 
%  species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
%  hold on
%  for i = 1:length(ts)
%      Sol = ys{i};
%      plot(ts{i}, sum(Sol(SumIndex1,:)),'r')
%      plot(ts{i}, sum(Sol(SumIndex2,:)),'g')
%      plot(ts{i}, sum(Sol(SumIndex3,:)),'b')
%      plot(ts{i}, sum(Sol(SumIndex4,:)),'c')
%      plot(ts{i}, sum(Sol(SumIndex5,:)),'m')
%      plot(ts{i}, sum(Sol(SumIndex6,:)),'k')
%  end
%  legend(species_in_graph); 

 
%  figure; 
%  species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo'};
%  
%  semilogx(Sol.x, sum(Sol.y(SumIndex1,:)))
%  hold on
%  semilogx(Sol.x, sum(Sol.y(SumIndex2,:)))
%  hold on
%  semilogx(Sol.x, sum(Sol.y(SumIndex3,:)))
%  hold on
%  semilogx(Sol.x, sum(Sol.y(SumIndex4,:)))
%  hold on 
%  semilogx(Sol.x, sum(Sol.y(SumIndex5,:)))
%  hold on
%  semilogx(Sol.x, sum(Sol.y(SumIndex6,:)))
% 
%  legend(species_in_graph); 
% 
% 
%  Fl = 1/(1+k(kn)+k(kr)+k(kq))*(Sol.y(y(YoPoAo),:)+Sol.y(y(YoPoAoBoo),:)...
%      +Sol.y(y(YoPoAoBro),:)+Sol.y(y(YoPoAoBrr),:)+Sol.y(y(YoPoAr),:)...
%      +Sol.y(y(YoPoArBoo),:)+Sol.y(y(YoPoArBro),:)+Sol.y(y(YoPoArBrr),:)...
%      +1/(1+k(kp)+k(kn)+k(kr))*(Sol.y(y(YoPrAo),:)+Sol.y(y(YoPrAoBoo),:)...
%      +Sol.y(y(YoPrAoBro),:)+Sol.y(y(YoPrAoBrr),:))+1/(1+k(kn)+k(kr))*(Sol.y(y(YoPrAr),:)...
%      +Sol.y(YoPrAoBoo)+Sol.y(y(YoPrArBro),:)+Sol.y(y(YoPrArBrr),:)...
%      +1/(1+k(kp)+k(kn))*(Sol.y(y(YrPrAo),:)+Sol.y(y(YrPrAoBoo),:)+Sol.y(y(YrPrAoBro),:)...
%      +Sol.y(YoPrArBrr))+1/(1+k(kn))*Sol.y(y(YrPrAr),:)+Sol.y(y(YrPrArBoo),:)...
%      +Sol.y(y(YrPrArBro),:)+Sol.y(y(YrPrArBrr),:)));

%{
 figure; 
 species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo'};
 
 semilogx(t, sum(S(SumIndex1,:)))
 hold on
 semilogx(t, sum(S(SumIndex2,:)))
 hold on
 semilogx(t, sum(S(SumIndex3,:)))
 hold on
 semilogx(t, sum(S(SumIndex4,:)))
 hold on 
 semilogx(t, sum(S(SumIndex5,:)))
 hold on
 semilogx(t, sum(S(SumIndex6,:)))

 legend(species_in_graph); 
 
 
 figure('Name','Y-Vals')
 species_in_graph = {'YoPrArBoo','YoPrAoBoo','YrPrAo'};
 idcs = [];
 
 for i = 1:length(species_in_graph)
     idcs(i) = find(strcmp(species,species_in_graph{i}));  
 end
 
 plot(Sol.x,Sol.y(idcs,:))
 legend(species_in_graph);

 figure('Name','Rates')
 YrPrAoBooidc = find(strcmp(species,'YrPrAoBoo'));
 indcs_of_reactions = find(S(YrPrAoBooidc,:)); 
 
 plot(Sol.x,r(indcs_of_reactions,:))
 Reactions = {};
 Reactions = (Rknames(indcs_of_reactions));
 legend(Reactions);
 
 
 
figure;
plot(Sol.x,sum(dydt));
legend('dydtSum');
foo = 1;
 
figure
species_in_graph = {'YrPrAo', 'YoPrAr', 'YrPrAr','YoPrAo', 'YoPoAr', 'YoPoAo' };
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(Sol.x,Sol.y(idcs,:))
legend(species_in_graph);

% YrPrAoBooidc = find(strcmp(species,'YrPrAoBoo'));
% indcs_of_reactions = find(S(YrPrAoBooidc,:)); 
% idcs = [];
% 
% for i = 1:length(indcs_of_reactions)
%     idcs = find(S(:,indcs_of_reactions(i)));
%     k = [i,:] 
% end 
% 
% figure;
% 
% for i = 1:length(indcs_of_reactions)


% plot(Sol.x, sum(Sol.y(idcs,:)));
 
% 
% for i = 1:length(species_in_graph)
%     figure; hold on;
%     plot(Sol.x, Sol.y(idcs(i),:));
%     plot(Sol.x, sum(Sol.y(idcs,:)));
%     legend(Ynames(i))
% 
% end
 
foo = 1;
 
%TEST PARAMETERS FROM TABLE 1
%
%make new folders for each condition they are testing in the paper
%
%run and save the same figures as they do apart from fluorescence. Also do
%O2
 
 %}
end
 
 
 
 
 
 
 
 
 
 
 
 
 