function [ts,ys, Fs, FvFm, species, O2] = main_FvFm_ode15s1(analysis_name,randomseed)


if nargin == 1
    randomseed = 'shuffle';
    
end
rng(randomseed);

rng(randomseed);
file1 = [analysis_name,'/FvFm_exp.csv'];
C = readtable(file1);
FvFm_exp = C.FqFm;

file1 = [analysis_name,'/experimental_parameters.xls'];
tbl1 = readtable(file1);
n_flashes = tbl1.n_flashes;
flash_duration = tbl1.flash_duration;
flash_interval = tbl1.flash_interval;
train_interval = tbl1.train_interval;
n_trains = tbl1.n_trains;


file1 = [analysis_name,'/LaiskConstants.xls'];
tablek = readtable(file1);
indepk = find(tablek.independent); 
lbk = tablek.lb(indepk);
ubk = tablek.ub(indepk);
knames = tablek.name;
k_std = tablek.base_val;
kwf = find(strcmp(knames,'kwf'));
kwb = find(strcmp(knames,'kwb'));

file2 = [analysis_name,'/LaiskY.xls'];
tabley = readtable(file2);
indepy = find(tabley.independent);
y_std = tabley.base_val;
lby = tabley.lb(indepy);
uby = tabley.ub(indepy);
yr = lby + (uby-lby).*rand(length(indepy),1);
kwfindep = length(yr) + length(find(tablek.independent(1:kwf)));
kwbindep = length(yr) + length(find(tablek.independent(1:kwb))); 
nspecies  = length(tabley.base_val);
deltapsiindex = nspecies + 1;
yidcs.deltapsiindex = deltapsiindex;
yinitial(deltapsiindex) = 0;
ATPaseoindex = nspecies + 2; 
yidcs.ATPaseoindex = ATPaseoindex;
yinitial(ATPaseoindex) = .63;
ATPaserindex = nspecies + 3; 
yidcs.ATPaserindex = ATPaserindex;
yinitial(ATPaserindex) = 0; 
pH_lumenindex = nspecies + 4;
yidcs.pH_lumenindex = pH_lumenindex;
yinitial(pH_lumenindex) = 10^-7.8;
pH_stromaindex = nspecies + 5;
yidcs.pH_stromaindex = pH_stromaindex;
yinitial(pH_stromaindex) = 10^-7.8;
fRindex = nspecies + 6;
yidcs.fRindex = fRindex;
yinitial(fRindex) = 0; 

% yr = yr+ yr.*rand(length(yr),1)*.5 - yr.*rand(length(yr),1)*.5;
kx = k_std(indepk);
yr = y_std(indepy);
% kx = lbk + (ubk-lbk).*rand(length(indepk),1);



indep_ps2_cell = cellfun(@(x) regexp(x, regexptranslate('wildcard', 'S*Y*P*')), tabley.name(indepy), 'UniformOutput', false);
indep_ps2 = [];
for i = 1:length(indep_ps2_cell)
    if ~isempty(indep_ps2_cell{i})
        indep_ps2(end+1) = i;
    end
end
Ay = zeros(1,length(yr));
Ay(indep_ps2) = 1;
by = 1;

Ak = zeros(1,length(indepk));
Aeq = [Ay,Ak];
beq = 1;

x0 = [yr;kx];

lb = [reshape(lby,[],1); reshape(lbk,[],1)];
ub = [reshape(uby,[],1); reshape(ubk,[],1)];

Ynames = tabley.name;
file3 = [analysis_name,'/LaiskReactions.xlsx'];
[~,Rknames] = xlsread(file3);

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
 
mult1 = find(strcmp(knames,'n2*kp/(1+kp+kn+kr)'));
mult2 = find(strcmp(knames,'n2*kp/(1+kp+kn)')); 
Div1 = find(strcmp(knames,'kpc/kEpc'));
Div2 = find(strcmp(knames,'kcytf/kEcytf'));
Div3 = find(strcmp(knames,'kfx/kEfx'));
Div4 = find(strcmp(knames,'kb6f/kEb6f'));
n1idx = find(strcmp(knames,'n1'));

kidcs.PFD = PFD;
kidcs.Labs = Labs;
kidcs.a2 = a2;
kidcs.PSU2 = PSU2;
kidcs.Chl = Chl;
kidcs.PSU1 = PSU1;
kidcs.kp = kp;
kidcs.kn = kn;
kidcs.kr = kr;
kidcs.Div1 = Div1;
kidcs.Div2 = Div2;
kidcs.Div3 = Div3;
kidcs.Div4 = Div4;
kidcs.mult1 = mult1;
kidcs.mult2 = mult2;
kidcs.n1idx = n1idx;
kidcs.kpc = kpc;
kidcs.kEpc = kEpc;
kidcs.kcytf = kcytf;
kidcs.kEcytf = kEcytf;
kidcs.kfx = kfx;
kidcs.kEfx = kEfx;
kidcs.kb6f = kb6f;
kidcs.kEb6f = kEb6f;
kidcs.kf = kf;

Sth = 1.56e2; 
F = 96485;
Cmem = 6e-3;
kpsi = 10;
T = 298; 
R = 8.31; 
apmf = 3e-2;
bpmf = 1.87; 


ATPpar.Thr = 0; 
ATPpar.Tho = .1; 
ATPpar.Stroma = 25e-6; 
ATPpar.EmATPase_7 = -.27; 
ATPpar.EmATPTh_7 = -.27;
ATPpar.kFr = 10e2;
ATPpar.B_stroma = 3.33e-2; 
ATPpar.Vlumen = 8.09e-7;
ATPpar.Vstroma = 6.48e-6; 
ATPpar.B_lumen = 3.33e-2;
ATPpar.pmfd = 6e-2; 
ATPpar.kF20 = 2.16e3; 
ATPpar.kF10 = 5.13e3;
ATPpar.kFC = 3.1;
ATPpar.VpH = 59; 
ATPpar.kF = 1.1e2; 
ATPpar.Sth = Sth;
ATPpar.F = F;
ATPpar.Cmem = Cmem;
ATPpar.kpsi = kpsi;
ATPpar.T = T;
ATPpar.R = R;
ATPpar.apmf = apmf;
ATPpar.bpmf = bpmf;

[species,S,rate_inds] = Laisk_read_excel_model1(analysis_name);

FDP = find(strcmp(species,'SynFDP'));
FP = find(strcmp(species,'SynFP'));
FD = find(strcmp(species,'SynFD'));
FT = find(strcmp(species,'SynFT'));
F = find(strcmp(species,'SynF'));
Hs = find(strcmp(species,'Hs'));
Hl = find(strcmp(species,'Hl'));
NADPH = find(strcmp(species,'NADPH'));
NADP = find(strcmp(species,'NADP'));
ADP = find(strcmp(species,'SynADP'));
ATP = find(strcmp(species,'SynATP'));
% Hs -> Hl 	kHleak
% ATP -> ADP 	kcbb_ATP
% NADPH -> NADP	kcbb_NADPH
% 9ATP + 6 NADPH -> 9 ADP + 6 NADP	kcbb

yidcs.FDP = FDP;
yidcs.FT = FT;
yidcs.FD = FD;
yidcs.FP = FP;
yidcs.F = F;
yidcs.Hs = Hs;
yidcs.Hl = Hl; 
yidcs.NADPH = NADPH;
yidcs.NADP = NADP;
yidcs.ADP = ADP;
yidcs.ATP = ATP; 

species_idcs = zeros(length(species),1);
% yinitial = zeros(length(y0),1);
for i = 1:length(Ynames)
    index = find(strcmp(species,Ynames(i)));
%     yinitial(index) = y0(i);
    species_idcs(index) = i;
end

PSIidcs = zeros(2,1);
PSIidcs(1) = find(strcmp(species,'P700o'));
PSIidcs(2) = find(strcmp(species,'FXo'));

[kconst] = LaiskKconstantsReadTable(analysis_name);

kn = find(strcmp(knames,'kn'));
kp = find(strcmp(knames,'kp'));
kr = find(strcmp(knames,'kr')); 
kq = find(strcmp(knames,'kq')); 
Fluorescence_k_idcs = [kn;kp;kr;kq];

kf1index = find(strcmp(knames,'kF1')); 
kf1indcs = find(kconst == kf1index);
kf2index = find(strcmp(knames,'kF2')); 
kf2indcs = find(kconst == kf2index); 

yopoax = find(contains(species,'YoPoA'));
yoprao = find(contains(species,'YoPrAo'));
yoprar = find(contains(species,'YoPrAr'));
yrprao = find(contains(species,'YrPrAo'));
yrprar = find(contains(species,'YrPrAr'));

Fluorescence_y_inds = {yopoax;yoprao;yoprar;yrprao;yrprar};

[ts,ys, Fs, FvFm, O2, O2_light, O2_dark] = calc_Species_concs_ode15s1(x0,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
                    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
                    Fluorescence_k_idcs, Fluorescence_y_inds,... indeces used to calculate fluorescence
                    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
                    tablek, tabley,... information on the k and y variables
                    kconst, rate_inds, S, species, knames, species_idcs, Rknames, analysis_name,yidcs,ATPpar,kf1indcs, kf2indcs); % model specific variables

save([analysis_name,'/result.mat'], 'ts','ys','Fs','FvFm','species','O2')
                
end





