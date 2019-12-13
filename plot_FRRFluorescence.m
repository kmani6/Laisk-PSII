function plot_FRRFluorescence(analysis_name)
load([analysis_name, '/FRR_results.mat'])
Print_Params(analysis_name)

PQ = find(strcmp(species, 'PQ'));
PQH2 = find(strcmp(species, 'PQH2'));
PCo = find(strcmp(species, 'PCo'));
PCr = find(strcmp(species, 'PCr'));
FDo = find(strcmp(species, 'FDo'));
FDr = find(strcmp(species, 'FDr'));
Hs = find(strcmp(species, 'Hs'));
Hl = find(strcmp(species, 'Hl'));
ATP = find(strcmp(species, 'ATP'));
ADP = find(strcmp(species, 'ADP'));
NADP = find(strcmp(species, 'NADP'));
NADPH = find(strcmp(species, 'NADPH'));


Fo = [];
Fm = [];
figure(1);
hold on
title('PQ redox state')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    plot(t, y(PQH2,:)./(y(PQH2,:) + y(PQ,:)),'b')
    
end

figure(2)
hold on
title('PC redox state')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    plot(t, y(PCr,:)./(y(PCr,:) + y(PCo,:)),'b')
end
figure(3)
hold on
title('Fd redox state')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    
    plot(t, y(FDr,:)./(y(FDr,:) + y(FDo,:)),'b')
end
figure(4)
hold on
title('Fluorescence')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    if ~isempty(Fs{i})
        F = Fs{i};
%         figure(4)
%         hold on
        plot(t,F,'b');
        Fo(end+1) = F(1);
        Fm(end+1) = F(end);        
    end
    %         figure(5)
    %         hold on
    %         plot(t(1:end-1),diff(F)./diff(t),'b')
end
FvFm = (Fm - Fo)./Fm;
figure(5)
plot(1:length(FvFm), FvFm,'-o')
title('FvFm')

figure(6)
hold on
plot(1:length(Fm), Fm, '-o');
plot(1:length(Fm), Fo, '-o');
legend({'Fm','Fo'})


figure(7)
hold on
title('ATP - ADP')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    plot(t, y(ATP,:),'b')
    plot(t, y(ADP,:),'k')
end
legend({'ATP','ADP'})

figure(8)
hold on
title('NADP - NADPH')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    plot(t, y(NADP,:),'b')
    plot(t, y(NADPH,:),'k')
end
legend({'NADP','NADPH'})

figure(9)
hold on
title('Hs - Hl')
for i = 2:length(ts)
    t = ts{i};
    y = ys{i};
    plot(t, y(Hs,:),'b')
    plot(t, y(Hl,:),'k')
end
legend({'Hs','Hl'})
end



