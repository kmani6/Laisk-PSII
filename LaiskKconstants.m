function [kconst] = LaiskKconstants(analysis_name)

file1 = [analysis_name,'/LaiskConstants.xls'];
[~,knames] = xlsread(file1);

file3 = [analysis_name,'/LaiskReactions.xls'];
[~,Rknames] = xlsread(file3);

jd = find(strcmp(knames,'jd')); 
oqd = find(strcmp(knames,'oqd'));
rqr = find(strcmp(knames,'rqr'));
rqd = find(strcmp(knames,'rqd'));
b1d = find(strcmp(knames,'b1d'));
b2d = find(strcmp(knames,'b2d'));
oqr = find(strcmp(knames,'oqr'));
b1r = find(strcmp(knames,'b1r'));
mult1 = find(strcmp(knames,'n2*kp/(1+kp+kn+kr)'));
mult2 = find(strcmp(knames,'n2*kp/(1+kp+kn)')); 
div1 = find(strcmp(knames,'kpc/kEpc'));
div2 = find(strcmp(knames,'kcytf/kEcytf'));
div3 = find(strcmp(knames,'kfx/kEfx'));
div4 = find(strcmp(knames,'kb6f/kEb6f'));
kcytf = find(strcmp(knames,'kcytf'));
kpc = find(strcmp(knames,'kpc'));
n1 = find(strcmp(knames, 'n1'));
kfx = find(strcmp(knames,'kfx'));
kfd = find(strcmp(knames,'kfd'));
kb6f = find(strcmp(knames,'kb6f'));

N = length(Rknames);
kconst = zeros(N,1);
 
rqd1 = find(strcmp(Rknames(:,2),'rqd'));
rqr1 = find(strcmp(Rknames(:,2),'rqr')); 
jd1 = find(strcmp(Rknames(:,2),'jd')); 
oqd1 = find(strcmp(Rknames(:,2),'oqd')); 
oqr1 = find(strcmp(Rknames(:,2),'oqr'));
b1d1 = find(strcmp(Rknames(:,2),'b1d')); 
b2d1 = find(strcmp(Rknames(:,2),'b2d')); 
b1r1 = find(strcmp(Rknames(:,2),'b1r')); 
Mult1 = find(strcmp(Rknames(:,2),'n2*kp/(1+kp+kn+kr)'));
Mult2 = find(strcmp(Rknames(:,2),'n2*kp/(1+kp+kn)')); 
Div1 = find(strcmp(Rknames(:,2),'kpc/kEpc'));
Div2 = find(strcmp(Rknames(:,2),'kcytf/kEcytf'));
Div3 = find(strcmp(Rknames(:,2),'kfx/kEfx'));
Div4 = find(strcmp(Rknames(:,2),'kb6f/kEb6f'));
kcytf1 = find(strcmp(Rknames(:,2),'kcytf')); 
kpc1 = find(strcmp(Rknames(:,2),'kpc'));
N1 = find(strcmp(Rknames(:,2),'n1'));
kfx1 = find(strcmp(Rknames(:,2),'kfx'));
kfd1 = find(strcmp(Rknames(:,2),'kfd'));
kb6f1 = find(strcmp(Rknames(:,2),'kb6f'));

% mult = find(all(isnan(x),2))

kconst(jd1) = jd;
kconst(oqd1) = oqd;
kconst(rqr1) = rqr;
kconst(rqd1) = rqd;
kconst(oqr1) = oqr;
kconst(b1d1) = b1d;
kconst(b2d1) = b2d;
kconst(b1r1) = b1r;
kconst(Mult1) = mult1;
kconst(Mult2) = mult2; 
kconst(Div1) = div1;
kconst(Div2) = div2;
kconst(Div3) = div3;
kconst(Div4) = div4; 
kconst(kcytf1) = kcytf;
kconst(kpc1) = kpc; 
kconst(N1) = n1;
kconst(kfx1) = kfx;
kconst(kfd1) = kfd;
kconst(kb6f1) = kb6f; 

end