function dydt = LaiskODEs(t,y,k,knames,Ynames)

PFD = find(strcmp(knames,'PFD'));
a2 = find(strcmp(knames,'a2'));
b1d = find(strcmp(knames,'b1d'));
b2d = find(strcmp(knames,'b2d'));
Chl = find(strcmp(knames,'Chl'));
CytfT = find(strcmp(knames,'CytfT'));
FDT = find(strcmp(knames,'FDT'));
jd = find(strcmp(knames,'jd')); 
kb6f = find(strcmp(knames,'kb6f'));
kcytf = find(strcmp(knames,'kcytf'));
kf = find(strcmp(knames,'kf'));
kfd = find(strcmp(knames,'kfd'));
kfx = find(strcmp(knames,'kfx'));
kn = find(strcmp(knames,'kn'));
kp = find(strcmp(knames,'kp'));
kpc = find(strcmp(knames,'kpc'));
kr = find(strcmp(knames,'kr'));
kE1 = find(strcmp(knames,'kE1'));
kEb6f = find(strcmp(knames,'kEb6f'));
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
         
PS1T = 1; 
PS2T = k(a2)*k(Chl)/k(PSU2);
n1 = k(PFD)*k(Labs)*(1-k(a2))/PS1T;
n2 = k(PFD)*k(Labs)*(k(a2))/PS2T; 
k(P700T) = PS1T; 
k(FXT) = PS1T; 

y(PCo) = k(PCT) - y(PCr); 
y(P700o) = k(P700T) - y(P700r); 
y(FXo) = k(FXT) - y(FXr); 
y(FDo) = k(FDT) - y(FDr);
y(Cytfo) = k(CytfT) - y(Cytfr);
Vb6f = k(kb6f)*(y(PQH2)*y(Cytfo)-y(Cytfr)*(y(PQ)/k(kEb6f)));

dydt = zeros(length(Ynames),1);
 
dydt(YoPoAo) = -k(jd)*y(YoPoAo)+k(oqd)*y(YoPoAoBoo) + k(rqd)*y(YoPoAoBrr) - (k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YoPoAo); %check
dydt(YoPoAoBoo) = -k(jd)*y(YoPoAoBoo)-k(oqd)*y(YoPoAoBoo) + k(oqr)*y(PQ)*y(YoPoAo);				%check	
dydt(YoPrAo) = k(jd)*y(YoPoAo) - k(jd)*y(YoPrAo) + k(rqd)* y(YoPrAoBrr) + k(oqd)*y(YoPrAoBoo) - (k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YoPrAo) - n2*k(kp)/(1+k(kp)+k(kn)+k(kr))*y(YoPrAo);	%check	
dydt(YoPoAr) = -k(jd)*y(YoPoAr) + n2*k(kp)/(1+k(kp)+k(kn)+k(kr)) * y(YoPrAo) + k(rqd)*y(YoPoArBrr)+k(oqd)*y(YoPoArBoo)-(k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YoPoAr);		%check
dydt(YoPrAoBoo) = -k(jd)*y(YoPrAoBoo) + k(jd)*y(YoPoAoBoo) - n2*(kp/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBoo)-k(oqd)*y(YoPrAoBoo)+k(oqr)*y(PQ)*y(YoPrAo);		%check
dydt(YoPoArBoo) = -k(jd)*y(YoPoArBoo)+n2*(k(kp)/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBoo)+k(oqr)*y(PQ)*y(YoPoAr)-k(oqd)*y(YoPoArBoo)-k(b1d)*y(YoPoArBoo)+k(b1r)*y(YoPoAoBro);		%check
dydt(YoPoAoBro) = -k(jd)*y(YoPoAoBro) + k(b1d)*y(YoPoArBoo)-k(b1r)*y(YoPoAoBro);						
dydt(YrPrAo) = k(jd)*y(YoPrAo)-n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAo)+k(rqd)*y(YrPrAoBrr)+k(oqd)*y(YrPrAoBoo)-(k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YrPrAo);		
dydt(YoPrAr) = k(jd)*y(YoPoAr)-k(jd)*y(YoPrAr)+n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAo)+k(rqd)*y(YoPrArBrr)+k(oqd)*y(YoPrArBoo)-(k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YoPrAr);		
dydt(YrPrAoBoo) = k(jd)*y(YoPrAoBoo)-n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBoo)+k(oqr)*y(PQ)*y(YrPrAo)-k(oqd)*y(YrPrAoBoo); %check		
dydt(YoPrArBoo) = k(jd)*y(YoPoArBoo)-k(jd)*y(YoPrArBoo)+n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBoo)-k(b1d)*y(YoPrArBoo)+k(b1r)*y(YoPrAoBro)-k(oqd)*y(YoPrArBoo)+k(oqr)*y(PQ)*y(YoPrAr);		
dydt(YoPrAoBro) = k(jd)*y(YoPoAoBro)-k(jd)*y(YoPrAoBro)-n2*(k(kp)/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBro)+k(b1d)*y(YoPrArBoo)-k(b1r)*y(YoPrAoBro);		
dydt(YoPoArBro) = -k(jd)*y(YoPoArBro)+n2*(k(kp)/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBro)-k(b2d)*y(YoPoArBro);		
dydt(YoPoAoBrr) = -k(jd)*y(YoPoAoBrr)+k(b2d)*y(YoPoArBro)+k(rqr)*y(PQH2)*y(YoPoAo)-k(rqd)*y(YoPoAoBrr);		 		
dydt(YrPrAr) = k(jd)*y(YoPrAr)+k(rqd)*y(YrPrArBrr)+k(oqd)*y(YrPrArBoo)-(k(oqr)*y(PQ)+k(rqr)*y(PQH2))*y(YrPrAr);		
dydt(YrPrArBoo) = k(jd)*y(YoPrArBoo)-k(b1d)*y(YrPrArBoo)+k(b1r)*y(YrPrAoBro)+k(oqr)*y(PQ)*y(YrPrAr)-k(oqd)*y(YrPrArBoo);		
dydt(YoPrArBro) = k(jd)*y(YoPoArBro)-k(jd)*y(YoPrArBro)+n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBro)-k(b2d)*y(YoPrArBro);		
dydt(YoPrAoBrr) = k(jd)*y(YoPoAoBrr)-k(jd)*y(YoPrAoBrr)-n2*(k(kp)/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBrr)+k(b2d)*y(YoPrArBro)-k(rqd)*y(YoPrAoBrr)+k(rqr)*y(PQH2)*y(YoPrAo);		
dydt(YrPrAoBro) = k(jd)*y(YoPrAoBro)-n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBro)+k(b1d)*y(YrPrArBoo)-k(b1r)*y(YrPrAoBro);		
dydt(YoPoArBrr) = -k(jd)*y(YoPoArBrr)+n2*(k(kp)/(1+k(kp)+k(kn)+k(kr)))*y(YoPrAoBrr)-k(rqd)*y(YoPoArBrr)+k(rqr)*y(PQH2)*y(YoPoAr);			
dydt(YrPrArBro) = k(jd)*y(YoPrArBro)-k(b2d)*y(YrPrArBro);		
dydt(YrPrAoBrr) = k(jd)*y(YoPrAoBrr)-n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBrr)+k(b2d)*y(YrPrArBro)-k(rqd)*y(YrPrAoBrr)+k(rqr)*y(PQH2)*y(YrPrAo);		
dydt(YoPrArBrr) = k(jd)*y(YoPoArBrr)-k(jd)*y(YoPrArBrr)+n2*(k(kp)/(1+k(kp)+k(kn)))*y(YrPrAoBrr)-k(rqd)*y(YoPrArBrr)+k(rqr)*y(PQH2)*y(YoPrAr);				
dydt(YrPrArBrr) = k(jd)*y(YoPrArBrr)-k(rqd)*y(YrPrArBrr)+k(rqr)*y(PQH2)*y(YrPrAr);		
dydt(PQH2) = k(rqd)*(y(YoPoAoBrr)+y(YoPrAoBrr)+y(YoPoArBrr)+y(YrPrAoBrr)+y(YoPrArBrr)+y(YrPrArBrr))-k(rqr)*y(PQH2)*(y(YoPoAo)+y(YoPrAo)+y(YoPoAr)+y(YrPrAo)+y(YoPrAr)+y(YrPrAr)) - Vb6f;
dydt(PQ) = k(oqd)*(y(YoPoAoBoo)+y(YoPrAoBoo)+y(YoPoArBoo)+y(YrPrAoBoo)+y(YoPrArBoo)+y(YrPrArBoo))-k(oqr)*y(PQ)*(y(YoPoAo)+y(YoPrAo)+y(YoPoAr)+y(YrPrAo)+y(YoPrAr)+y(YrPrAr))+ Vb6f;
dydt(Cytfr) = 2*Vb6f-k(kcytf)*(y(Cytfr)*y(PCo)-y(Cytfo)*(y(PCr)/k(kEcytf)));
dydt(PCr) = k(kcytf)*(y(Cytfr)*y(PCo)-y(Cytfo)*(y(PCr)/k(kEcytf)))-k(kpc)*(y(PCr)*y(P700o)-y(PCo)*(y(P700r)/k(kEpc)));
dydt(P700r) = k(kpc)*(y(PCr)*y(P700o)-y(PCo)*(y(P700r)/kEpc)) - n1*y(P700r)*y(FXo);
dydt(FXr) = n1*y(P700r)*y(FXo)-k(kfx)*(y(FXr)*y(FDo)-y(FXo)*(y(FDr)/k(kEfx)));
dydt(FDr) = k(kfx)*(y(FXr)*y(FDo)-y(FXo)*(y(FDr)/k(kEfx))) - k(kfd)*y(FDr);
dydt(O2) = k(jd)*(y(YoPoAo)+y(YoPoAoBoo)+y(YoPrAo)+y(YoPoAr)+y(YoPrAoBoo)+y(YoPoArBoo)+y(YoPoAoBro)+y(YoPrAr)+y(YoPrArBoo)+y(YoPrAoBro)+y(YoPoArBro)+y(YoPoAoBrr)+y(YoPrArBro)+y(YoPrAoBrr)+y(YoPoArBrr)+y(YoPrArBrr));

Fl = (1/(1+k(kn)+k(kr)+k(kq)))*((y(YoPoAo)+y(YoPoAoBoo)+y(YoPoAoBro)+y(YoPoAoBrr)+y(YoPoAr)+y(YoPoArBoo)+y(YoPoArBro)+y(YoPoArBrr))+(1/(1+k(kp)+k(kn)+k(kr)))*((y(YoPrAo)+y(YoPrAoBoo)+y(YoPrAoBro)+y(YoPrAoBrr)))+(1/(1+k(kn)+k(kr)))*((y(YoPrAr)+y(YoPrArBoo)+y(YoPrArBro)+y(YoPrArBrr))) + (1/(1+k(kp)+k(kn)))*((y(YrPrAo)+y(YrPrAoBoo)+y(YrPrAoBro)+y(YrPrArBrr)))+(1/(1+k(kn)))*((y(YrPrAr)+y(YrPrArBoo)+y(YrPrArBro)+y(YrPrArBrr))));

dydt(F) = n2*(1-Fl);


if t>0
    foo = 1;
end
end

%dydt(PCo) = -dydt(PCr); 
%dydt(P700o) = -dydt(P700r); 
%dydt(FXo) = -dydt(FXr);
%dydt(FDo) = -dydt(FDr);
%dydt(Cytfo) = -dydt(Cytfr); 



