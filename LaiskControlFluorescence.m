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
n1 = k(PFD)*k(Labs)*(1-k(a2))/PS1T;
n2 = k(PFD)*k(Labs)*k(a2)/PS2T; 
y0 = y0*PS2T; 

PS1T = 1; 
PS2T = 1;

k(P700T) = PS1T; 
k(FXT) = PS1T; 

% y(PCo) = k(PCT) - y(PCr); 
% y(P700o) = k(P700T) - y(P700r); 
% y(FXo) = k(FXT) - y(FXr); 
% y(FDo) = k(FDT) - y(FDr);
% y(Cytfo) = k(CytfT) - y(Cytfr);

%Vb6f = k(kb6f)*(y(PQH2)*y(Cytfo)-y(Cytfr)*(y(PQ)/k(kEb6f)));

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

[kconst] = LaiskKconstants(analysis_name);

Sol = ode15s(@(t,y) LaiskPS2ODES(t,y,kconst,rate_inds,S),[tstart,tend],y0);
dydt = [];
for i = 1:length(Sol.x)
    dydt(:,i) = LaiskPS2ODES(Sol.x(i),Sol.y(:,i),kconst,rate_inds,S);
end

foo = 1;

%[Fl] = fluorescence(Ynames,knames,k,Sol);

%figure(1)
%semilogx(Sol.x,Fl)

%hold on 
%plot(Sol.x,Sol.y(end,:))
%plot(Sol.x,Sol.y(end-1,:))


figure;
species_in_graph = Ynames(1:24);
idcs = [];

for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(species,species_in_graph{i}));
   
end
semilogx(Sol.x,Sol.y(idcs,:))

figure;
semilogx(Sol.x, sum(Sol.y(idcs,:)))
    title('YxPxAxBxx');
legend(species_in_graph)
figure
species_in_graph = {'PQ','PQH2'};
idcs = [];
for i = 1:length(Ynames)  
    idcs(i) = find(strcmp(Ynames,species_in_graph{i}));

end 
    semilogx(Sol.y,Sol.y(idcs,:))
    legend(Ynames(i))
figure
species_in_graph = {'PQ','PQH2'};
idcs = [];
hold on
for i = 1:length(species_in_graph)
    idcs(i) = find(strcmp(Ynames,species_in_graph{i}));
    semilogx(Sol.x,Sol.y(idcs(i),:))
end

semilogx(Sol.x, sum(Sol.y(idcs,:)))
legend('PQ','PQH2');
title('PQ');


Cytfr = find(strcmp(Ynames,'Cytfr'));
CytfT = find(strcmp(knames,'CytfT'));
figure
hold on

plot(Sol.x, Sol.y(Cytfr,:))
plot(Sol.x, k(CytfT)*PS2T-Sol.y(Cytfr,:))
legend('Ctyfr','Cytfo');
title('Cytf');

hold off

PCr = find(strcmp(Ynames,'PCr'));
PCT = find(strcmp(knames,'PCT'));
figure
hold on

plot(Sol.x, Sol.y(PCr,:))
plot(Sol.x, k(PCT)-Sol.y(PCr,:))
legend('PCr','PCo');
title('PC');

hold off

P700r = find(strcmp(Ynames,'P700r'));
a2 = find(strcmp(knames,'a2'));
Chl = find(strcmp(knames,'Chl'));
PSU1 = find(strcmp(knames,'PSU1'));
PS1T = 1; %(1-k(a2))*(k(Chl)/k(PSU1));
P700T = find(strcmp(knames, 'P700T'));
FXT = find(strcmp(knames, 'FXT')); 
k(P700T) = PS1T; 
k(FXT) = PS1T; 
FDT = find(strcmp(knames,'FDT'));

figure
hold on

plot(Sol.x, Sol.y(P700r,:))
plot(Sol.x, k(P700T)-Sol.y(P700r,:))
legend('P700r','P700o');
title('P700');

hold off

FXr = find(strcmp(Ynames,'FXr')); 

figure
hold on

plot(Sol.x, Sol.y(FXr,:))
plot(Sol.x, k(FXT)-Sol.y(FXr,:))
legend('FXr','FXo');
title('FX');

hold off

FDr = find(strcmp(Ynames,'FDr')); 

figure
hold on

plot(Sol.x, Sol.y(FDr,:))
plot(Sol.x, k(FDT)-Sol.y(FDr,:))
legend('FDr','FDo');
title('FD');

hold off

end 

% figure 
% semilogx(Sol.x, Sol.y(YoPoAo,:))
% ylim([-10^-6,10^-6])
% title('YoPoAo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoAoBoo,:))
% ylim([-10^-6,10^-6])
% title('YoPoAoBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrAo,:))
% ylim([-10^-6,10^-6])
% title('YoPrAo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoAr,:))
% ylim([-10^-6,10^-6])
% title('YoPoAr');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrAoBoo,:))
% ylim([-10^-6,10^-6])
% title('YoPrAoBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoArBoo,:))
% ylim([-10^-6,10^-6])
% title('YoPoArBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoAoBro,:))
% ylim([-10^-6,10^-6])
% title('YoPoAoBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrAo,:))
% ylim([-10^-6,10^-6])
% title('YrPrAo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrAr,:))
% ylim([-10^-6,10^-6])
% title('YoPrAr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrAoBoo,:))
% ylim([-10^-6,10^-6])
% title('YrPrAoBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrArBoo,:))
% ylim([-10^-6,10^-6])
% title('YoPrArBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrAoBro,:))
% ylim([-10^-6,10^-6])
% title('YoPrAoBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoArBro,:))
% ylim([-10^-6,10^-6])
% title('YoPoArBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoAoBrr,:))
% ylim([-10^-6,10^-6])
% title('YoPoAoBrr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrAr,:))
% ylim([-10^-6,10^-6])
% title('YrPrAr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrArBoo,:))
% ylim([-10^-6,10^-6])
% title('YrPrArBoo');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrArBro,:))
% ylim([-10^-6,10^-6])
% title('YoPrArBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrAoBrr,:))
% ylim([-10^-6,10^-6])
% title('YoPrAoBrr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrAoBro,:))
% ylim([-10^-6,10^-6])
% title('YrPrAoBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPoArBrr,:))
% ylim([-10^-6,10^-6])
% title('YoPoArBrr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrArBro,:))
% ylim([-10^-6,10^-6])
% title('YrPrArBro');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrAoBrr,:))
% ylim([-10^-6,10^-6])
% title('YrPrAoBrr');
% 
% figure
% semilogx(Sol.x, Sol.y(YoPrArBrr,:))
% ylim([-10^-6,10^-6])
% title('YoPrArBrr');
% 
% figure
% semilogx(Sol.x, Sol.y(YrPrArBrr,:))
% ylim([-10^-6,10^-6])
% title('YrPrArBrr');