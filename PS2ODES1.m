function dydt = PS2ODES1(t,y,krxn,k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs)
% disp(num2str(t))

nrxn = length(rate_inds);
fFr = (y(yidcs.ATPaserindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
fFo = (y(yidcs.ATPaseoindex))/(y(yidcs.ATPaserindex) + y(yidcs.ATPaseoindex));
deltaph = y(yidcs.pH_lumenindex) - y(yidcs.pH_stromaindex);
pmf = y(yidcs.deltapsiindex) + (2.3*ATPpar.R*ATPpar.T)*(deltaph)/ATPpar.F; 
x = (1/ATPpar.kF)*(fFr*10^pmf/ATPpar.VpH + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH); 
D = 1 + x + x^2 + x^3 + x^4; 
p1 = (ATPpar.kFC*x^4)/((1+ATPpar.kFC)*D); 
p2 = 1/((1+ATPpar.kFC)*D); 
krxn(kf1indcs) = ATPpar.kF10*p1; 
krxn(kf2indcs) = ATPpar.kF20*p2; 

r = zeros(nrxn,1);
if any(isnan(y))
    foo = 1;
end
for irxn = 1:nrxn
    r(irxn,1) = krxn(irxn)*prod(y(rate_inds{irxn}));
end

%deltaph = -log(y(yidcs.Hs)) + log(y(yidcs.Hl)); 
r(268,1) = r(268,1)*(fFr*10^(pmf/ATPpar.VpH) + fFo*10^(pmf-ATPpar.pmfd)/ATPpar.VpH);

dydt = S*r;
% dFl = 0;% dLaiskFluorescence(Ynames,knames,k,y);
% dydt(end+1) = dFl;
pmf_ss = min([ATPpar.apmf + (ATPpar.bpmf*2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph, (ATPpar.apmf+ATPpar.bpmf)*2.3*ATPpar.R*ATPpar.T/ATPpar.F * deltaph]);
deltapsi_ss = (pmf_ss)*(-2.3*ATPpar.R*ATPpar.T/ATPpar.F)*deltaph; 
kpsi = 10;
Fi = (deltapsi_ss-y(yidcs.deltapsiindex))*kpsi*ATPpar.Cmem/ATPpar.F;
VHin = dydt(yidcs.Hl);
VHout = 0; %dydt(yidcs.Hs); 
dydt(yidcs.deltapsiindex) = (((VHin-VHout)/ATPpar.Sth + Fi)*ATPpar.F)/ATPpar.Cmem;
%Ftotal = y(yidcs.FDP) + y(yidcs.FP) + y(yidcs.FDP) + y(yidcs.FT); 

Em_Th = ATPpar.EmATPTh_7 *ATPpar.VpH*(-log(y(yidcs.Hs)-7));
Em_ATPase = ATPpar.EmATPase_7*ATPpar.VpH*(-log(y(yidcs.Hs)-7));
KE_ThATPase = exp((2*ATPpar.F*(Em_ATPase-Em_Th))/(ATPpar.R*ATPpar.T));
VFr = (ATPpar.kFr/ATPpar.Stroma)*(ATPpar.Thr*y(yidcs.ATPaseoindex)-((ATPpar.Tho*y(yidcs.ATPaserindex))/(KE_ThATPase)));
dydt(yidcs.ATPaserindex) = ATPpar.B_stroma*VFr;
dydt(yidcs.ATPaseoindex) = -ATPpar.B_stroma*VFr;
dydt(yidcs.pH_stromaindex) = (ATPpar.B_stroma/ATPpar.Vstroma) - dydt(yidcs.Hs);
dydt(yidcs.pH_lumenindex) = (ATPpar.B_lumen/ATPpar.Vlumen) - dydt(yidcs.Hl); 

%fFo = y(yidcs.FT)/Ftotal; 
%fFo initial val = 0 



end


