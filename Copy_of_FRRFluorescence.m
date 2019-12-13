function FvFm = FRRFluorescence(analysis_name)

file1 = [analysis_name,'/experimental_parameters'];
[parameters, parameter_names] = xlsread(file1);
n_flashes = parameters(1);
flash_duration = parameters(2);
flash_interval = parameters(3);
train_interval = parameters(4);
n_trains = parameters(5);

FvFm = zeros(n_flashes*n_trains,1);

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
dark_adaptation_time = 5; %3-5 minutes typically
t = linspace(0, dark_adaptation_time, dark_adaptation_time*1e5);
t_lims = [0,dark_adaptation_time];
tic;
Sol =  ode15s(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t_lims,yinitial);
toc
% Sol = Sol';
ys{end+1} = Sol.y; %calculate the species evolution during the light
ts{end+1} = -dark_adaptation_time + Sol.x; %save the times used
Fs{end+1} = []; %save an empty fluorescence vector
yinitial = Sol.y(:,end); %initialize the y vector for the next iteration 

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
counter = 1
for train = 1:n_trains
    fprintf('train %i \n', train)
    for flash = 1:n_flashes
        k(mult1) = mult1Val;
        k(mult2) = mult2Val;    
        k(n1idx) = n1;
        t = linspace(0, flash_duration, flash_duration*1e6);
        Sol = ode2(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
        Sol = Sol';
        ys{end+1} = Sol; %calculate the species evolution during the light
        ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + (flash-1)*(flash_duration+flash_interval)+t; %save the times used
        F = LaiskFluorescence(species,knames,k,Sol); 
        Fs{end+1} = F; %save the fluoprescence values
        yinitial = Sol(:,end); %initialize the y vector for the next iteration 
        FvFm(counter) = (max(F) - F(1))/max(F);
        counter = counter+1;
        %Shift to dark time between flashes
        k(mult1) = 0;
        k(mult2) = 0;      
        k(n1idx) = 0;
        t = linspace(0, flash_interval, flash_interval*1e6);
        t_lims = [0,flash_interval];
        Sol = ode15s(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t_lims,yinitial);
        ys{end+1} = Sol.y; %calculate the species evolution during the dark between flashes
        ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + (flash-1)*(flash_duration+flash_interval) + flash_duration +Sol.x; %save the times used
        Fs{end+1} = []; % save an empty vector for fluorescence (to allign the values of times and fluorescence)
        yinitial = Sol.y(:,end); %initialize the y vector for the next iteration 
        
    end
    k(mult1) = 0;
    k(mult2) = 0;       
    k(n1idx) = 0;
    t = linspace(0, train_interval, train_interval*1e5);
    Sol = ode15s(@(t,y) PS2ODES(y,k(kconst),k,rate_inds,S,species,knames),t,yinitial);
    ys{end+1} = Sol.y; %calculate the species evolution during the dark between flashes
    ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+Sol.x; %save the times used
    Fs{end+1} = []; % save an empty vector for fluorescence (to allign the values of times and fluorescence)
    yinitial = Sol.y(:,end); %initialize the y vector for the next iteration    
end

save([analysis_name '/FRR_results.mat'], 'ts', 'Fs', 'ys', 'knames', 'kconst', 'k', 'species', 'Rknames','-v7.3')  

 
%TEST PARAMETERS FROM TABLE 1
%
%make new folders for each condition they are testing in the paper
%
%run and save the same figures as they do apart from fluorescence. Also do
%O2
 
 %}
end
 
 
 
 
 
 
 
 
 
 
 
 
 