function dydt = PS2ODES1(t,y,krxn,k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs, kidcs)
% disp(num2str(t))

nrxn = length(rate_inds);
fFr = (y(yidcs.ATPaserindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
fFo = 1-fFr;%(y(yidcs.ATPaseoindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
deltaph = y(yidcs.pH_stromaindex) - y(yidcs.pH_lumenindex); %y(yidcs.pH_lumenindex) - y(yidcs.pH_stromaindex);
deltapsi = .59*deltaph;
y(yidcs.deltapsiindex) = deltapsi;
pmf = y(yidcs.deltapsiindex)+(2.3*ATPpar.R*ATPpar.T)*(deltaph)/ATPpar.F; 
x = (1/ATPpar.kF)*10^(-y(yidcs.pH_lumenindex))/10^(-y(yidcs.pH_stromaindex));
% x = (1/ATPpar.kF)*(fFr*10^(pmf/ATPpar.VpH) + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH); 
D = 1 + x + x^2 + x^3 + x^4; 
p1 = (ATPpar.kFC*x^4)/((1+ATPpar.kFC)*D); 
p2 = 1/((1+ATPpar.kFC)*D); 
krxn(kf1indcs) = ATPpar.kF10*p1*100; 
krxn(kf2indcs) = ATPpar.kF20*p2; 
r = zeros(nrxn,1);
if t>0
    foo = 1;
end
for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end
dydt = S*r;

if t > 1e-8
   foo = 1;  
end 


HPR = 4.67; 
phi = 1; 
vFDPFT = (y(yidcs.FDP)*krxn(kf1indcs)-y(yidcs.FT)*krxn(kf2indcs)); 
Vrmax = 12;
KmPGA = 3.5; 
KmNADPH = 0.023;
KmATP = 0.24;
PGA = 19;
fRmin = 0.04;
fRalpha = 0.00248; %in m^2*s/umol
fRtheta = 0.96;
fRss = fRmin + ((fRalpha * k(kidcs.PFD) + (1 - fRmin)) - sqrt((fRalpha * k(kidcs.PFD) + (1 - fRmin))^2 - 4 * fRtheta * fRalpha * k(kidcs.PFD) * (1 - fRmin))) / (2 * fRtheta);
kiR = 6.3e-3;
kdR = 7.5e-3; 
if fRss > y(yidcs.fRindex)
        dydt(yidcs.fRindex) = (fRss-y(yidcs.fRindex))*kiR;
else
        dydt(yidcs.fRindex) = (fRss-y(yidcs.fRindex))*kdR;
end 
Vr = (PGA*y(yidcs.ATP)*y(yidcs.NADPH)*y(yidcs.fRindex)*Vrmax)/((PGA+KmPGA)*(y(yidcs.ATP)+KmATP)*(y(yidcs.NADPH)+KmNADPH)); 
dydt(yidcs.ATP) = dydt(yidcs.ATP) - (3+3.5*phi)*Vr/(2+1.5*phi);
dydt(yidcs.ADP) = dydt(yidcs.ADP) + (3+3.5*phi)*Vr/(2+1.5*phi);
dydt(yidcs.NADPH) = dydt(yidcs.NADPH) - (2+2*phi)*Vr/(2+1.5*phi);
dydt(yidcs.NADP) = -dydt(yidcs.NADPH);
% a = dydt(yidcs.Hl);
aa = dydt(yidcs.Hl)/ATPpar.Vlumen;
aaaVfdpdt =vFDPFT;
aaa = vFDPFT * HPR/3/ATPpar.Vlumen;
aaat = t;
dydt(yidcs.Hl) = aa-aaa; 

% dydt(yidcs.Hl) = dydt(yidcs.Hl)/ATPpar.Vstroma - vFDPFT/ATPpar.Vstroma * HPR/3/ATPpar.Vstroma; 
dydt(yidcs.Hs) = dydt(yidcs.Hs)/ATPpar.Vstroma - ((2.0 + 2.0 * phi) * Vr) / (2 + 1.5 * phi) +(HPR * vFDPFT) / ATPpar.Vstroma;
% dydt(yidcs.Hl) = dydt(yidcs.Hl) - vFDPFT * HPR;% / ATPpar.Vlumen; 

% disp(a);
% dFl = 0;% dLaiskFluorescence(Ynames,knames,k,y);
% dydt(end+1) = dFl;
pmf_ss = min([ATPpar.apmf + (ATPpar.bpmf*2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph, (ATPpar.apmf+ATPpar.bpmf)*2.3*ATPpar.R*ATPpar.T/ATPpar.F * deltaph]);
deltapsi_ss = (pmf_ss)*(2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph; 
kpsi = 10;
v_counter = (deltapsi_ss - y(yidcs.deltapsiindex)) * kpsi;
Fi = (deltapsi_ss-y(yidcs.deltapsiindex))*kpsi*ATPpar.Cmem/ATPpar.F;
VHin = dydt(yidcs.Hl);
VHout = 0; %dydt(yidcs.Hs); 



dydt(yidcs.deltapsiindex) = dydt(yidcs.Hl)* ATPpar.Vlumen/ATPpar.thylac*ATPpar.F/ATPpar.Cmem+ v_counter;%(((VHin-VHout)/ATPpar.Sth + Fi)*ATPpar.F)/ATPpar.Cmem; %*1e-4;
%Ftotal = y(yidcs.FDP) + y(yidcs.FP) + y(yidcs.FDP) + y(yidcs.FT); 

Em_Th = ATPpar.EmATPTh_7-ATPpar.VpH*(y(yidcs.pH_stromaindex)-7);
Em_ATPase = ATPpar.EmATPase_7-ATPpar.VpH*(y(yidcs.pH_stromaindex)-7);
KE_ThATPase = exp((2*ATPpar.F*(Em_ATPase-Em_Th))/(ATPpar.R*ATPpar.T));
VFr = (ATPpar.kFr/ATPpar.Stroma)*(ATPpar.Thr*y(yidcs.ATPaseoindex)-((ATPpar.Tho*y(yidcs.ATPaserindex))/(KE_ThATPase)));
% dydt(yidcs.ATPaserindex) = ATPpar.B_stroma*VFr;
% dydt(yidcs.ATPaseoindex) = -ATPpar.B_stroma*VFr;
dydt(yidcs.pH_stromaindex) = -(ATPpar.B_stroma/ATPpar.Vstroma)*dydt(yidcs.Hs);
dydt(yidcs.pH_lumenindex) = -(ATPpar.B_lumen/ATPpar.Vlumen)*dydt(yidcs.Hl); 
% dydt(yidcs.FT) = dydt(yidcs.FT) - (r(kf1indcs)-r(kf2indcs)) + (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FDP) = dydt(yidcs.FDP) + (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FT) = dydt(yidcs.FT) - (r(kf1indcs)-r(kf2indcs)) + (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);
% dydt(yidcs.FDP) = dydt(yidcs.FDP) + (r(kf1indcs)-r(kf2indcs)) - (r(kf1indcs)-r(kf2indcs))*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);

%fFo = y(yidcs.FT)/Ftotal; 
%fFo initial val = 0 

if pmf>5
    foo = 1;
end

if t > 1e-8
   foo = 1;  
end 

if any(y<-1e-3)
   badies = find(y<-1e-3);
   y(badies)
   species(badies)
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
end


end


