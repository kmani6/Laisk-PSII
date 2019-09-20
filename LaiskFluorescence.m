function[Fl] = LaiskFluorescence(Ynames,knames,k,Sol) 

kn = find(strcmp(knames,'kn'));
kp = find(strcmp(knames,'kp'));
kr = find(strcmp(knames,'kr')); 
kq = find(strcmp(knames,'kq')); 

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
YoPrArBrr = find(strcmp(Ynames,'YoPrArBrr')); 
YrPrArBrr = find(strcmp(Ynames,'YrPrArBrr')); 
YoPrArBoo = find(strcmp(Ynames, 'YoPrArBoo'));
YrPrAoBrr = find(strcmp(Ynames, 'YrPrAoBrr'));

Fl = 1/(1+k(kn)+k(kr)+k(kq))*Sol(YoPoAo,:)+Sol(YoPoAoBoo,:)...
     +Sol(YoPoAoBro,:)+Sol(YoPoAoBrr,:)+Sol(YoPoAr,:)...
     +Sol(YoPoArBoo,:)+Sol(YoPoArBro,:)+Sol(YoPoArBrr,:)...
     +1/(1+k(kp)+k(kn)+k(kr))*(Sol(YoPrAo,:)+Sol(YoPrAoBoo,:)...
     +Sol(YoPrAoBro,:)+Sol(YoPrAoBrr,:))+1/(1+k(kn)+k(kr))*(Sol(YoPrAr,:)...
     +Sol(YoPrArBoo,:)+Sol(YoPrArBro,:)+Sol(YoPrArBrr,:)...
     +1/(1+k(kp)+k(kn))*(Sol(YrPrAo,:)+Sol(YrPrAoBoo,:)...
     +Sol(YrPrAoBro,:)+Sol(YrPrAoBrr))+1/(1+k(kn))*Sol(YrPrAr,:)...
     +Sol(YrPrArBoo,:)+Sol(YrPrArBro,:)+Sol(YrPrArBrr,:));

end
