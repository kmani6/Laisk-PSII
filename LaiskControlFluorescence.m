function[Sol] = LaiskControlFluorescence(tspan,analysis_name)
 
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
CytfT = find(strcmp(knames, 'Chl')); 
FDT = find(strcmp(knames, 'Chl')); 
jd = find(strcmp(knames, 'Chl')); 
kb6f = find(strcmp(knames, 'Chl')); 
kcytf = find(strcmp(knames, 'Chl')); 
kf = find(strcmp(knames, 'Chl')); 
kfd = find(strcmp(knames, 'Chl')); 
kfx = find(strcmp(knames, 'Chl')); 
kn = find(strcmp(knames, 'Chl')); 
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
kq = find(strcmp(knames, 'kq')); 
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
F = find(strcmp(Ynames, 'F')); 
 
PS1T = k(a2)*k(Chl)/k(PSU2); 
PS2T = (1-k(a2))*(k(Chl)/k(PSU1));
% PS1T = 1; 
% PS2T = 1;
n1 = k(PFD)*k(Labs)*(1-k(a2))/PS1T;
n2 = k(PFD)*k(Labs)*k(a2)/PS2T; 
y0 = y0*PS2T; 
 
 
k(P700T) = PS1T; 
k(FXT) = PS1T; 
 
Vb6f = k(kb6f)*(y0(PQH2)*y0(Cytfo)-y0(Cytfr)*(y0(PQ)/k(kEb6f)));
 
mult1 = find(strcmp(knames,'n2*kp/(1+kp+kn+kr)'));
mult2 = find(strcmp(knames,'n2*kp/(1+kp+kn)')); 
Div1 = find(strcmp(Rknames(:,2),'kpc/kEpc'));
Div2 = find(strcmp(Rknames(:,2),'kcytf/kEcytf'));
Div3 = find(strcmp(Rknames(:,2),'kfx/kEfx'));
Div4 = find(strcmp(Rknames(:,2),'kb6f/kEb6f'));
 
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
 
tstart = tspan(1);
tend = tspan(2);
 
[species,S,rate_inds] = Laisk_read_excel_model('Laisk DCMU.1');
yinitial = zeros(length(y0),1);
for i = 1:length(Ynames)
    index = find(strcmp(species,Ynames(i)));
    yinitial(index) = y0(i);
end
 
[kconst] = LaiskKconstants(analysis_name);
 
Sol = ode23t(@(t,y) LaiskPS2ODES(t,y,kconst,rate_inds,S),linspace(tstart,tend,1e3),yinitial);
dydt = [];
for i = 1:length(Sol.x)
    dydt(:,i) = LaiskPS2ODES(Sol.x(i),Sol.y(:,i),kconst,rate_inds,S);
end
  
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

semilogx(Sol.x,Sol.y(idcs,:))
legend(species_in_graph);

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
 
 
end
 
 
 
 
 
 
 
 
 
 
 
 
 