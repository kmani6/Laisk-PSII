function plot_PS2_species(Sol, t, species);

yopoax = find(contains(species,'YoPoA'));
yoprao = find(contains(species,'YoPrAo'));
yoprar = find(contains(species,'YoPrAr'));
yrprao = find(contains(species,'YrPrAo'));
yrprar = find(contains(species,'YrPrAr'));

figure;
hold on
plot(t, Sol(yopoax,:))
plot(t, Sol(yoprao,:))
plot(t, Sol(yoprar,:))
plot(t, Sol(yrprao,:))
plot(t, Sol(yrprar,:))

legend({'YoPoAx','YoPrAo','YoPrAr','YrPrAo','YrPrAr'})