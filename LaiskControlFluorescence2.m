function[Sol] = LaiskControlFluorescence2(tspan,analysis_name)
 
maxtime = tspan(end);
tmax  = tspan(2); 
 
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



tstart = tspan(1);
tend = tspan(2);
 
[species,S,rate_inds] = Laisk_read_excel_model(analysis_name);
yinitial = zeros(length(y0),1);
for i = 1:length(Ynames)
    index = find(strcmp(species,Ynames(i)));
    yinitial(index) = y0(i);
end

[kconst] = LaiskKconstants(analysis_name);
 
Fo = LaiskFluorescence(species,knames,k,yinitial);
ytest = yinitial;
ytest(YrPrAoBoo) = 0;
ytest(YrPrArBrr) = PS2T;
Fm =  LaiskFluorescence(species,knames,k,ytest);
% Sol = ode15s(@(t,y) LaiskPS2ODES(t,y,k(kconst),rate_inds,S),[tstart,tend],yinitial);
%Sol.x(1) = Sol.x(2)/100;
t = linspace(0, tend, 50000);
% t = logspace(-5, log10(tend), 5000);
t(1) = 0;
% oqr_inds = find(kconst == oqr);
% rqr_inds = find(kconst == rqr);
% PQ_idx = find(strcmp(species,'PQ'));
% PQH2_idx = find(strcmp(species,'PQH2'));
PQ = find(strcmp(species, 'PQ'));
PQH2 = find(strcmp(species, 'PQH2'));
rqr1 = find(strcmp(Rknames(:,2),'rqr')); 
oqr1 = find(strcmp(Rknames(:,2),'oqr'));
Sol = ode2(@(t,y) LaiskPS2ODES(t,y,k(kconst),k,rate_inds,S,species,knames,PQ,PQH2,oqr1,rqr1),t,yinitial);
Sol = Sol';
dydt = [];
for i = 1:length(t)
%     dydt(:,i) = LaiskPS2ODES(t,Sol(:,i),k(kconst),k,rate_inds,S,species,knames);
    r(:,i) = LaiskRates(t,Sol(:,i),k(kconst),rate_inds,S);
end


% for i =length(Rknames):-1:1
%     figure;
%     plot(t,r(i,:))
%     legend(Rknames(i))
%     
% end
% 
% 


 SumIndex1 = find(contains(species, 'YrPrAo'));
 SumIndex2 = find(contains(species, 'YoPrAr'));
 SumIndex3 = find(contains(species, 'YrPrAr'));
 SumIndex4 = find(contains(species, 'YoPrAo'));
 SumIndex5 = find(contains(species, 'YoPoAr'));
 SumIndex6 = find(contains(species, 'YoPoAo'));
 
 YoPoArBrr = find(contains(species, 'YoPoArBrr')); 
 figure;
 plot(t, Sol(YoPoArBrr,:),'o','color','red','MarkerSize',1.5);
 
 
 
 
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
Fl = LaiskFluorescence(species,knames,k,Sol); 
 
 
%  figure; 
%  species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
%  
%  semilogx(t, sum(Sol(SumIndex1,:)))
%  hold on
%  semilogx(t, sum(Sol(SumIndex2,:)))
%  hold on
%  semilogx(t, sum(Sol(SumIndex3,:)))
%  hold on
%  semilogx(t, sum(Sol(SumIndex4,:)))
%  hold on 
%  semilogx(t, sum(Sol(SumIndex5,:)))
%  hold on
%  semilogx(t, sum(Sol(SumIndex6,:)))
%  hold on
%  semilogx(t, Fl);
% %  hold on
% %  semilogx(t, Sol(end,:));
%  
%  legend(species_in_graph); 
 

%  figure; 
%  species_in_graph = {'YrPrAo','YoPrAr','YrPrAr','YoPrAo','YoPoAr','YoPoAo','Fl'};
%  
%  plot(t, sum(Sol(SumIndex1,:)))
%  hold on
%  plot(t, sum(Sol(SumIndex2,:)))
%  hold on
%  plot(t, sum(Sol(SumIndex3,:)))
%  hold on
%  plot(t, sum(Sol(SumIndex4,:)))
%  hold on 
%  plot(t, sum(Sol(SumIndex5,:)))
%  hold on
%  plot(t, sum(Sol(SumIndex6,:)))
%  hold on
%  plot(t, Fl);
% %  hold on
% %  semilogx(t, Sol(end,:));
%  
%  legend(species_in_graph); 
 

% figure;
% hold on
% species_in_graph = {{'PQH2','PQ'}, {'Cytfr','Cytfo'} {'PCr','PCo'}, {'P700r','P700o'}, {'FDr','FDo'}};
% lgd = {};
% for i = 1:length(species_in_graph)
%     
%     current_species = species_in_graph{i};
%     idcs = [];
%     for j = 1:length(current_species)
%         idcs(end+1) = find(strcmp(species, current_species{j}));        
%     end
%     redox_state = Sol(idcs(1),:)./(Sol(idcs(1),:) + Sol(idcs(2),:));
%     plot(t,redox_state,'o','MarkerSize',3.5);
%     lgd{end+1} = (strcat(current_species{1}, "/(", current_species{2}, " + ", current_species{1}, ")"))  ;
% end
% plot(t,Fl,'o','MarkerSize',3.5);
% set(gca,'FontSize',20)
% set(gca,'linewidth',2)
% ylim([0 1.2])
% lgd{end+1} = "Fl";
% legend(lgd);
figure; 
i1 = find(strcmp(species,'P700o'));
i2 = find(strcmp(species,'P700r'));
y1 = Sol(i1,:);
y2 = Sol(i2,:);
plot(t, [y1; y2]); 
legend('P700o', 'P700r')
 foo = 1;


figure;
species_in_graph = {'PQ', 'PQH2'};
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(t,Sol(idcs,:))
legend(species_in_graph);


figure;
species_in_graph = {'PCr', 'PCo'};
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(t,Sol(idcs,:))
legend(species_in_graph);

figure;
species_in_graph = {'Cytfr', 'Cytfo'};
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(t,Sol(idcs,:))
legend(species_in_graph);

figure;
species_in_graph = {'FXr', 'FXo'};
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(t,Sol(idcs,:))
legend(species_in_graph);

figure;
species_in_graph = {'FDr', 'FDo'};
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));  
end

plot(t,Sol(idcs,:))
legend(species_in_graph);

foo = 1;
%  
% 
%  
% figure;
% species_in_graph = species(find(contains(species,'Y')));
% idcs = [];
% 
% for i = 1:length(species_in_graph)
%     idcs(i) = find(strcmp(species,species_in_graph{i}));  
% end
% 
% plot(t,Sol(idcs,:))
% legend(species_in_graph);
%  
% rs = [];



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

save([analysis_name,'/result.mat'], 'Sol','t','species','k','knames','kconst','Ynames','y0','Rknames')

end
 
 
 