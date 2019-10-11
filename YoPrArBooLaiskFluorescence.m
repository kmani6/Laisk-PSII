function[Fl1] = YoPrArBooLaiskFluorescence(Ynames,knames,k,Sol) 

kn = find(strcmp(knames,'kn'));
kr = find(strcmp(knames,'kr')); 

YoPrArBoo = find(strcmp(Ynames,'YoPrArBoo')); 

Fl1 = 1/(1+k(kn)+k(kr))*(Sol(YoPrArBoo,:));
    
end
