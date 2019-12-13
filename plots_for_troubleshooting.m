y = ys{2};
t = ts{2};
f = Fs{2};
pq = find(strcmp(species,'PQ'));
pqh2 = find(strcmp(species,'PQH2'));
figure; hold on; plot(t,y(pqh2,:)./(y(pq,:) + y(pqh2,:)))
pco = find(strcmp(species,'PCo'));
pcr = find(strcmp(species,'PCr'));
plot(t,y(pcr,:)./(y(pco,:) + y(pcr,:)))
io = find(strcmp(species,'FDo'));
ir = find(strcmp(species,'FDr'));
plot(t,y(ir,:)./(y(io,:) + y(ir,:)))
io = find(strcmp(species,'NADP'));
ir = find(strcmp(species,'NADPH'));
plot(t,y(ir,:)./(y(io,:) + y(ir,:)))
plot(t,f)

legend({'PQ', 'PC', 'FD', 'NADP', 'Fl'})


figure; plot(t,y(pqh2,:),t,y(pq,:)); 
legend({'PQH2', 'PQ'})

boo = find(contains(species,'Boo'));
bro = find(contains(species,'Bro'));
brr = find(contains(species,'Brr'));
nob = find(~contains(species,'B') & contains(species, 'Y'));


figure; plot(t,sum(y(boo,:)))
hold on
plot(t,sum(y(bro,:)))
plot(t,sum(y(brr,:)))
plot(t,sum(y(nob,:)))
legend({'Boo','Bro','Brr','Empty site'})



yopoax = find(contains(species,'YoPoA'));
yoprao = find(contains(species,'YoPrAo'));
yoprar = find(contains(species,'YoPrAr'));
yrprao = find(contains(species,'YrPrAo'));
yrprar = find(contains(species,'YrPrAr'));


figure; plot(t,sum(y(yopoax,:)))
hold on
plot(t,sum(y(yoprao,:)))
plot(t,sum(y(yoprar,:)))
plot(t,sum(y(yrprao,:)))
plot(t,sum(y(yrprar,:)))
legend({'YoPoAx','YoPrAo','YoPrAr','YrPrAo','YrPrAr'})





s0 = find(contains(species,'S0'));
s1 = find(contains(species,'S1'));
s2 = find(contains(species,'S2'));
s3 = find(contains(species,'S3'));


figure; plot(t,sum(y(s0,:)))
hold on
plot(t,sum(y(s1,:)))
plot(t,sum(y(s2,:)))
plot(t,sum(y(s3,:)))
legend({'S0','S1','S2','S3'})