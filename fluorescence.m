function[Fl] = fluorescence(Ynames,knames,k,Sol)

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

Fl = 1/(1+k(kn)+k(kr)+k(kq))*(Sol(y(YoPoAo),:)+Sol(y(YoPoAoBoo),:)+Sol(y(YoPoAoBro),:)+Sol(y(YoPoAoBrr),:)+Sol(y(YoPoAr),:)+Sol(y(YoPoArBoo),:)+Sol(y(YoPoArBro),:)+Sol(y(YoPoArBrr),:)+1/(1+k(kp)+k(kn)+k(kr))*(Sol(y(YoPrAo),:)+Sol(y(YoPrAoBoo),:)+Sol(y(YoPrAoBro),:)+Sol(y(YoPrAoBrr),:))+1/(1+k(kn)+k(kr))*(Sol(y(YoPrAr),:)+Sol(YoPrAoBoo)+Sol(y(YoPrArBro),:)+Sol(y(YoPrArBrr),:) + 1/(1+k(kp)+k(kn))*(Sol(y(YrPrAo),:)+Sol(y(YrPrAoBoo),:)+Sol(y(YrPrAoBro),:)+Sol(YoPrArBrr))+1/(1+k(kn))*Sol(y(YrPrAr),:)+Sol(y(YrPrArBoo),:)+Sol(y(YrPrArBro),:)+Sol(y(YrPrArBrr),:)));

end
