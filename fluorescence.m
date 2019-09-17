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

Fl = 1/(1+k(kn)+k(kr)+k(kq))*(Sol.y(y(YoPoAo),:)+Sol.y(y(YoPoAoBoo),:)+Sol.y(y(YoPoAoBro),:)+Sol.y(y(YoPoAoBrr),:)+Sol.y(y(YoPoAr),:)+Sol.y(y(YoPoArBoo),:)+Sol.y(y(YoPoArBro),:)+Sol.y(y(YoPoArBrr),:)+1/(1+k(kp)+k(kn)+k(kr))*(Sol.y(y(YoPrAo),:)+Sol.y(y(YoPrAoBoo),:)+Sol.y(y(YoPrAoBro),:)+Sol.y(y(YoPrAoBrr),:))+1/(1+k(kn)+k(kr))*(Sol.y(y(YoPrAr),:)+Sol.y(YoPrAoBoo)+Sol.y(y(YoPrArBro),:)+Sol.y(y(YoPrArBrr),:) + 1/(1+k(kp)+k(kn))*(Sol.y(y(YrPrAo),:)+Sol.y(y(YrPrAoBoo),:)+Sol.y(y(YrPrAoBro),:)+Sol.y(YoPrArBrr))+1/(1+k(kn))*Sol.y(y(YrPrAr),:)+Sol.y(y(YrPrArBoo),:)+Sol.y(y(YrPrArBro),:)+Sol.y(y(YrPrArBrr),:)));

end
