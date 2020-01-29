function plot_Fl(Fs, ts)
figure;
% hold on;
t = [];
F = [];
for i = 1:length(ts)
    if ~isempty(Fs{i})
        t = [t;reshape(ts{i},[],1)];
        F = [F;reshape(Fs{i},[],1)];
    end
end
plot(t,F);
xlabel('time')
ylabel('Fl')
title('Fluorescence')

end