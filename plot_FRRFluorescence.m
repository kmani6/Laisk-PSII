function plot_FRRFluorescence(analysis_name)
load([analysis_name, '/FRR_results.mat'])

PQ = find(strcmp(species, 'PQ'));
PQH2 = find(strcmp(species, 'PQH2'));
PCo = find(strcmp(species, 'PCo'));
PCr = find(strcmp(species, 'PCr'));
FDo = find(strcmp(species, 'FDo'));
FDr = find(strcmp(species, 'FDr'));


figure(1)
title("PQ redox state")


figure(2)
title("PC redox state")

figure(3)
title("Fd redox state")

figure(4)
title("Fluorescence")

figure(5)
title("time")
counter = 1;
for i = 2:length(ts)
    t = ts{i};
    figure(5)
    hold on
    plot(counter:counter+length(t)-1, t,'b')
    counter = counter+length(t);
%     
%     y = ys{i};
%     
%     figure(1)
%     redox = y(PQH2,:)./(y(PQH2,:) + y(PQ,:));
%     hold on
%     plot(t, redox,'b')
%     figure(2)
%     redox = y(PCr,:)./(y(PCr,:) + y(PCo,:));
%     hold on
%     plot(t, redox,'b')
%     figure(3)
%     hold on
%     redox = y(FDr,:)./(y(FDr,:) + y(FDo,:));
%     plot(t, redox,'b')
%     clear redox
%     if ~isempty(Fs{i})
%         F = Fs{i};
%         figure(4)
%         hold on
%         plot(t,F,'b');
%     end
%     
end



