function dydt = PS2ODES1(t,y,krxn,k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs, kidcs)
% disp(num2str(t))
global PS2T
nrxn = length(rate_inds);

% fFr = (y(yidcs.ATPaserindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
% fFo = (1-fFr); %(y(yidcs.ATPaseoindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
% deltaph = y(yidcs.pH_stromaindex) - y(yidcs.pH_lumenindex); %y(yidcs.pH_lumenindex) - y(yidcs.pH_stromaindex);
% pmf = y(yidcs.deltapsiindex) + (2.3*ATPpar.R*ATPpar.T)*(deltaph)/ATPpar.F; 
% x = (1/ATPpar.kF)*(fFr*10^(pmf/ATPpar.VpH) + fFo*10^((pmf-ATPpar.pmfd)/ATPpar.VpH)); 
% D = 1 + x + x^2 + x^3 + x^4; 
% p1 = (ATPpar.kFC*x^4)/((1+ATPpar.kFC)*D); 
% p2 = 1/((1+ATPpar.kFC)*D); 
% krxn(kf1indcs) = ATPpar.kF10*p1; 
% krxn(kf2indcs) = ATPpar.kF20*p2; 
% 
% r = zeros(nrxn,1);
% if t>2.222e-05
%     foo = 1;
% end

PIt = 30;% in mM
Pi = PIt - 3*y(yidcs.ATP) - 2* y(yidcs.ADP) - 2.0*y(yidcs.RuBP) - y(yidcs.PGA);
y(yidcs.P) = Pi;

one_molar = 1000;
y(yidcs.Hs) = (10^-y(yidcs.pH_stromaindex)) * one_molar;
y(yidcs.Hl) = (10^-y(yidcs.pH_lumenindex)) * one_molar;

for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end

KFc = 3.1;
KF = 110;
kF10 = 5130.0;
kFT20 = 2160.0;
pmfd = 60;
V_per_pH = 59;
fFr = 1;
%fFo = 1-fFr;

dpHC = ((y(yidcs.pH_stromaindex)- y(yidcs.pH_lumenindex)) * 2.3 * ATPpar.R * ATPpar.T) / ATPpar.F;
pmf = y(yidcs.deltapsiindex) + dpHC;
x = (fFr * (10.0^(pmf / V_per_pH)) + (1 - fFr) * (10.0^((pmf - pmfd) / V_per_pH))) / KF;
D = 1.0 + x + (x^2.0) + (x^3.0) + (x^4.0);
p1 = ((KFc / (1.0 + KFc)) * (x^4.0)) / D;
p2 = (1.0 / (1.0 + KFc)) / D;
kF1 = kF10 * p1;
kFT2 = kFT20 * p2;
vFDPFT = y(yidcs.FDP) * kF1 - y(yidcs.FT) * kFT2;

r(265) = y(yidcs.FDP) * kF1;
r(266) = y(yidcs.FT) * kFT2;

dydt = S*r;

%(y(yidcs.ATPaserindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
%(y(yidcs.ATPaseoindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
% deltaph = -y(yidcs.pH_lumenindex) + y(yidcs.pH_stromaindex);
% fFo = (y(yidcs.ATPaseoindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
% deltaph = y(yidcs.pH_stromaindex) - y(yidcs.pH_lumenindex); %y(yidcs.pH_lumenindex) - y(yidcs.pH_stromaindex);
% y(yidcs.deltapsiindex) = 59*deltaph;
% pmf = y(yidcs.deltapsiindex) + (2.3*ATPpar.R*ATPpar.T)*(deltaph)/ATPpar.F; 
% x = (1/ATPpar.kF)*10^(-y(yidcs.pH_lumenindex))/(10^(-y(yidcs.pH_stromaindex)));
% x = (1/ATPpar.kF)*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH); 
% D = 1 + x + x^2 + x^3 + x^4; 
% p1 = (ATPpar.kFC*x^4)/((1+ATPpar.kFC)*D); 
% p2 = 1/((1+ATPpar.kFC)*D); 
% krxn(kf1indcs) = ATPpar.kF10*p1; 
% krxn(kf2indcs) = ATPpar.kF20*p2; 

HPR = 4.67;
lumen  = 1.30;
stroma  = 10.44;
thylakoid  = 251;

% PGA = 19;% in mM initial value
KmRuBP = 0.02;% in mM
KiPGA = 0.84;% in mM
Kmapp_RuBP = KmRuBP*(1 + y(yidcs.PGA)/KiPGA);% in mM

%RuBP = 1; in mM initial value

fRB = 0.4; %initial value;
kc = 4.4; % in 1/s
Kmc = 309; % in umol/mol
Kmo = 179; % in mmol/mol
RB = 1.45; % in mM # 1.6. Adjust to new volume
O2 = 210; % in mmol/mol
Cc = .4; % in mmol/mol
ko = 1; % in 1/s
fRuBP = (1 / (2 * RB)) * ((RB + Kmapp_RuBP + y(yidcs.RuBP)) - sqrt(((RB + Kmapp_RuBP + y(yidcs.RuBP))^2) - 4 * RB * y(yidcs.RuBP)));
Vc = (fRuBP * fRB * kc * RB * Cc) / (Cc + Kmc * (1.0 + O2 / Kmo));
Vo = (fRuBP * fRB * ko * RB * O2) / (O2 + Kmo * (1.0 + Cc / Kmc));
if Vc>0
    phi = Vo/Vc;
else
    phi = 0;
end

%Moved Before Dydt calculations
% x = (fFr * (10.0^(pmf / V_per_pH)) + (1 - fFr) * (10.0^((pmf - pmfd) / V_per_pH))) / KF;
% D = 1.0 + x + (x^2.0) + (x^3.0) + (x^4.0);
% p1 = ((KFc / (1.0 + KFc)) * (x^4.0)) / D;
% p2 = (1.0 / (1.0 + KFc)) / D;
% kF1 = kF10 * p1;
% kFT2 = kFT20 * p2;
% vFDPFT = y(yidcs.FDP) * kF1 - y(yidcs.FT) * kFT2;
% 
% r(265) = y(yidcs.FDP) * kF1;
% r(266) = y(yidcs.FT) * kFT2;
%------------------------------

%deltaph = -log(y(yidcs.Hs)) + log(y(yidcs.Hl)); 
% r(264,1) = r(264,1)*(fFr*10^(pmf/ATPpar.VpH) + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% r(265,1) = r(265,1)*(fFr*10^(pmf/ATPpar.VpH) + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);

% dydt = S*r;
% HPR = 4.67;
% phi = 1; 
% vFDPFT = (y(yidcs.FDP)*krxn(kf1indcs)-y(yidcs.FT)*krxn(kf2indcs)); 
% Vrmax = 12;
% KmPGA = 3.5; 
% KmNADPH = 0.023;
% KmATP = 0.24;
% PGA = 19;

fRmin = 0.04;
fRalpha = 0.00248;% in m^2*s/umol
fRtheta = 0.96;
fRss = fRmin + (fRalpha*k(kidcs.PFD) + (1 - fRmin) - sqrt( (fRalpha*k(kidcs.PFD) + (1 - fRmin))^2.0 - 4*fRtheta*fRalpha*k(kidcs.PFD)*(1 - fRmin) ))/(2*fRtheta);

fR_k_i = 6.3e-3;% in 1/s
fR_k_d = 7.5e-3;% in 1/s
fR = 0;

if fRss > fR
    dydt(yidcs.fRindex) = (fRss - y(yidcs.fRindex))*fR_k_i;
else
    dydt(yidcs.fRindex) =(fRss - y(yidcs.fRindex))*fR_k_d;
end

Vrmax = 12;% in mM/s # 10. Adjust new volume
KmPGA = 3.5;% in mM
KmNADPH = 0.023;% in mM
KmATP = 0.24;% in mM

Vr = y(yidcs.fRindex)*Vrmax*y(yidcs.PGA)/(y(yidcs.PGA) + KmPGA)*y(yidcs.NADPH)/(y(yidcs.NADPH) + KmNADPH)*y(yidcs.ATP)/(y(yidcs.ATP) + KmATP);% in mM/s

dydt(yidcs.PGA) = 2.0*Vc + 1.5*Vc*phi - Vr;
dydt(yidcs.RuBP) = (1.0 + phi)/(2.0 + 1.5*phi)*Vr - Vc*(1.0 + phi);% in mM/s


d_Hlt_dt = (-HPR * vFDPFT) / lumen ;
d_Hst_dt = (HPR * vFDPFT) / stroma - ((2.0 + 2.0 * phi) * Vr) / (2 + 1.5 * phi);

dydt(yidcs.Hl) = dydt(yidcs.Hl)/lumen + d_Hlt_dt; 
dydt(yidcs.Hs) = dydt(yidcs.Hs)/stroma + d_Hst_dt; 

Buffer_lumen = .0333;% in dm^3/mol
dydt(yidcs.pH_lumenindex) = -dydt(yidcs.Hl)*Buffer_lumen;% in 1/s


Buffer_stroma = .0333;% in dm^3/mol
dydt(yidcs.pH_stromaindex) = -dydt(yidcs.Hs)*Buffer_stroma;% in 1/s

pmf_a = 30;% in mV
pmf_b = 1.866;% in mV
pmfss = min([(dpHC * (pmf_a + pmf_b)),pmf_a + pmf_b * dpHC]);
DPsiSS = pmfss - dpHC;
K_PsiSS = 10;% in 1/s
v_counter = (DPsiSS - y(yidcs.deltapsiindex)) * K_PsiSS;
F = 96485.0;% in s*A/mol
Cmem = 0.6;%  in uF/cm^2
dydt(yidcs.deltapsiindex) = (((dydt(yidcs.Hl) * lumen) / thylakoid) / F) / Cmem + v_counter;


dydt(yidcs.ATP) = dydt(yidcs.ATP)  -(3 + 3.5*phi)*Vr/(2 + 1.5*phi);
dydt(yidcs.ADP) = dydt(yidcs.ADP)  + (3 + 3.5*phi)*Vr/(2 + 1.5*phi);
dydt(yidcs.NADPH) = dydt(yidcs.NADPH)  -(2.0 + 2.0*phi)*Vr/(2 + 1.5*phi);
dydt(yidcs.NADP) = -dydt(yidcs.NADPH);


if t > 1e-8
   foo = 1;  
end 


% HPR = 4.67; 
% phi = 1; 
% vFDPFT = (y(yidcs.FDP)*krxn(kf1indcs)-y(yidcs.FT)*krxn(kf2indcs)); 
% Vrmax = 12;
% KmPGA = 3.5; 
% KmNADPH = 0.023;
% KmATP = 0.24;
% PGA = 19;
% fRmin = 0.04;
% fRalpha = 0.00248; %in m^2*s/umol
% fRtheta = 0.96;
% fRss = fRmin + ((fRalpha * k(kidcs.PFD) + (1 - fRmin)) - sqrt((fRalpha * k(kidcs.PFD) + (1 - fRmin))^2 - 4 * fRtheta * fRalpha * k(kidcs.PFD) * (1 - fRmin))) / (2 * fRtheta);
% kiR = 6.3e-3;
% kdR = 7.5e-3; 
% if fRss > y(yidcs.fRindex)
%         dydt(yidcs.fRindex) = (fRss-y(yidcs.fRindex))*kiR;
% else
%         dydt(yidcs.fRindex) = (fRss-y(yidcs.fRindex))*kdR;
% end 
% Vr = (PGA*y(yidcs.ATP)*y(yidcs.NADPH)*y(yidcs.fRindex)*Vrmax)/((PGA+KmPGA)*(y(yidcs.ATP)+KmATP)*(y(yidcs.NADPH)+KmNADPH)); 
% dydt(yidcs.ATP) = dydt(yidcs.ATP)  -(3 + 3.5*phi)*Vr/(2 + 1.5*phi);
% dydt(yidcs.ADP) = dydt(yidcs.ADP)  + (3 + 3.5*phi)*Vr/(2 + 1.5*phi);
% dydt(yidcs.NADPH) = dydt(yidcs.NADPH)  -(2.0 + 2.0*phi)*Vr/(2 + 1.5*phi);
% dydt(yidcs.NADP) = -dydt(yidcs.NADPH);
% dydt(yidcs.Hs) = dydt(yidcs.Hs)/ATPpar.Vlumen + (HPR * vFDPFT)/ATPpar.Vlumen - ((2.0 + 2.0 * phi) * Vr) / (2 + 1.5 * phi);
% a = dydt(yidcs.Hl);
% aa = dydt(yidcs.Hl)/ATPpar.Vstroma;
% aaaVfdpdt =vFDPFT;
% aaa = vFDPFT/ATPpar.Vstroma * HPR/3;
% aaat = t;


% dydt(yidcs.Hl) = dydt(yidcs.Hl)/ATPpar.Vstroma - vFDPFT/ATPpar.Vstroma * HPR/3/ATPpar.Vstroma; 
% dydt(yidcs.Hs) = dydt(yidcs.Hs) + (-(HPR * vFDPFT) / ATPpar.Vstroma + ((2.0 + 2.0 * phi) * Vr) / (2 + 1.5 * phi));
% dydt(yidcs.Hl) = dydt(yidcs.Hl) + vFDPFT * HPR / ATPpar.Vlumen; 
% dydt(yidcs.Hs) = dydt(yidcs.Hs) + (-(HPR * vFDPFT) / ATPpar.Vstroma + ((2.0 + 2.0 * phi) * Vr) / (2 + 1.5 * phi));
% dydt(yidcs.Hl) = dydt(yidcs.Hl) + vFDPFT * HPR / ATPpar.Vlumen; 

% dFl = 0;% dLaiskFluorescence(Ynames,knames,k,y);
% dydt(end+1) = dFl;
% pmf_ss = min([ATPpar.apmf + (ATPpar.bpmf*2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph, (ATPpar.apmf+ATPpar.bpmf)*2.3*ATPpar.R*ATPpar.T/ATPpar.F * deltaph]);
% deltapsi_ss = (pmf_ss)*(2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph; 
% kpsi = 10;
% Fi = (deltapsi_ss-y(yidcs.deltapsiindex))*kpsi*ATPpar.Cmem/ATPpar.F;
% VHin = dydt(yidcs.Hl);
% VHout = 0; %dydt(yidcs.Hs); 
% dydt(yidcs.deltapsiindex) = (((VHin-VHout)/ATPpar.Sth + Fi)*ATPpar.F)/ATPpar.Cmem; %*1e-4;
% %Ftotal = y(yidcs.FDP) + y(yidcs.FP) + y(yidcs.FDP) + y(yidcs.FT); 

% Em_Th = ATPpar.EmATPTh_7-ATPpar.VpH*(y(yidcs.pH_stromaindex)-7);
% Em_ATPase = ATPpar.EmATPase_7-ATPpar.VpH*(y(yidcs.pH_stromaindex)-7);
% KE_ThATPase = exp((2*ATPpar.F*(Em_ATPase-Em_Th))/(ATPpar.R*ATPpar.T));
% VFr = (ATPpar.kFr/ATPpar.Stroma)*(ATPpar.Thr*y(yidcs.ATPaseoindex)-((ATPpar.Tho*y(yidcs.ATPaserindex))/(KE_ThATPase)));
% dydt(yidcs.ATPaserindex) = ATPpar.B_stroma*VFr;
% dydt(yidcs.ATPaseoindex) = -ATPpar.B_stroma*VFr;
% dydt(yidcs.pH_stromaindex) = -(ATPpar.B_stroma/ATPpar.Vstroma)*dydt(yidcs.Hs);
% dydt(yidcs.pH_lumenindex) = -(ATPpar.B_lumen/ATPpar.Vlumen)*dydt(yidcs.Hl); 
% dydt(yidcs.FDP) = dydt(yidcs.FDP) - (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.pH_stromaindex) = -(ATPpar.B_stroma)*dydt(yidcs.Hs);
% dydt(yidcs.pH_lumenindex) = -(ATPpar.B_lumen)*dydt(yidcs.Hl); 
% dydt(yidcs.FT) = dydt(yidcs.FT);% - (r(kf1indcs)-r(kf2indcs)) + (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FDP) = dydt(yidcs.FDP);% + (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FT) = dydt(yidcs.FT)- (r(kf1indcs)-r(kf2indcs)) + (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FDP) = dydt(yidcs.FDP) + (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% pmf_ss = min([ATPpar.apmf + (ATPpar.bpmf*2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph, (ATPpar.apmf+ATPpar.bpmf)*2.3*ATPpar.R*ATPpar.T/ATPpar.F * deltaph]);
% deltapsi_ss = (pmf_ss)*(-2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph; 
% kpsi = 10;
% Fi = (deltapsi_ss-y(yidcs.deltapsiindex))*kpsi*ATPpar.Cmem/ATPpar.F;
% VHin = dydt(yidcs.Hl);
% VHout = 0; %dydt(yidcs.Hs); 
% dydt(yidcs.deltapsiindex) = (((VHin-VHout)/ATPpar.Sth + Fi)*ATPpar.F)/ATPpar.Cmem; %*1e-4;
%Ftotal = y(yidcs.FDP) + y(yidcs.FP) + y(yidcs.FDP) + y(yidcs.FT); 

% Em_Th = ATPpar.EmATPTh_7 - ATPpar.VpH*(-log(y(yidcs.Hs)-7)); %Em_Th = ATPpar.EmATPTh_7 * ATPpar.VpH*(-log(y(yidcs.Hs)-7));
% Em_ATPase = ATPpar.EmATPase_7-ATPpar.VpH*(-log(y(yidcs.Hs)-7)); %ATPpar.EmATPase_7*ATPpar.VpH*(-log(y(yidcs.Hs)-7));
% KE_ThATPase = exp((2*ATPpar.F*(Em_ATPase-Em_Th))/(ATPpar.R*ATPpar.T));
% VFr = (ATPpar.kFr/ATPpar.Stroma)*(ATPpar.Thr*y(yidcs.ATPaseoindex)-((ATPpar.Tho*y(yidcs.ATPaserindex))/(KE_ThATPase)));
% dydt(yidcs.ATPaserindex) = ATPpar.B_stroma*VFr;
% dydt(yidcs.ATPaseoindex) = -ATPpar.B_stroma*VFr;
% <<<<<<< HEAD
% dydt(yidcs.pH_stromaindex) = -(ATPpar.B_stroma/ATPpar.Vstroma)*dydt(yidcs.Hs);
% dydt(yidcs.pH_lumenindex) = -(ATPpar.B_lumen/ATPpar.Vlumen)*dydt(yidcs.Hl); 
% dydt(yidcs.FDP) = dydt(yidcs.FDP) - (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% =======
% dydt(yidcs.pH_stromaindex) = -(ATPpar.B_stroma)*dydt(yidcs.Hs);
% dydt(yidcs.pH_lumenindex) = -(ATPpar.B_lumen)*dydt(yidcs.Hl); 
% >>>>>>> dfd6a071ed8da6ea4a02edd291f6cd46525d60c1
% dydt(yidcs.FT) = dydt(yidcs.FT) - (r(kf1indcs)-r(kf2indcs)) + (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FDP) = dydt(yidcs.FDP) + (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);

%fFo = y(yidcs.FT)/Ftotal; 
%fFo initial val = 0 
if t > .0009
   foo = 1;
end 

if any(y<0)
   badies = find(y<0);
   y(badies)
   species(badies)
%    Untitled2
   foo = 1;
end

if any(isnan(dydt))
   badies = find(isnan(dydt));
   dydt(badies)
   species(badies)
   foo = 1; 
end 
if any(isnan(y))
   badies = find(isnan(y));
   y(badies)
   species(badies)
   foo = 1; 
end 
if any(y>1e3)
   badies = find(y>1e3);
   y(badies)
   species(badies)
   foo = 1; 
% <<<<<<< HEAD
if any (y<-1e-10)
   foo = 1;
% =======
if any(y<-1e-10)
    foo = 1;
% >>>>>>> dfd6a071ed8da6ea4a02edd291f6cd46525d60c1
end
if isnan(dydt)
   foo = 1;
end
if isnan(y)
   foo = 1;
end

end 
end


