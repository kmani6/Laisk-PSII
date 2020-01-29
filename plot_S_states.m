function f = plot_S_states(species, ys, ts) % model specific variables)

s0 = contains(species,'S0');
s1 = contains(species,'S1');
s2 = contains(species,'S2');
s3 = contains(species,'S3');

f = figure;
t = [];
S0 = [];
S1 = [];
S2 = [];
S3 = [];
for time_per = 2:length(ys)
    S0 = [S0, sum(ys{time_per}(s0,:))];
    S1 = [S1, sum(ys{time_per}(s1,:))];
    S2 = [S2, sum(ys{time_per}(s2,:))];
    S3 = [S3, sum(ys{time_per}(s3,:))];
    t = [t, ts{time_per}];

end
hold on
plot(t,S0,'b')
plot(t,S1,'r')
plot(t,S2,'k')
plot(t,S3,'m')
legend({'S0', 'S1', 'S2', 'S3'})
title('S-state distribution')
ylabel('S-state fraction')
xlabel('time, s')
    
    
    

