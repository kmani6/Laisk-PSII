function [Sol,Fpulses,Fout] = LaiskFluorescence(tspan)

maxtime = tspan(end);
y0 = zeros(1,36);
k = zeros(1,32);
tmax = tspan(2);

%[kvals,knames] = xlsread(file);

file1 = ['Laisk','/LaiskConstants.xls'];
[kvals,knames] = xlsread(file1);

file2 = ['Laisk','/LaiskY.xls'];
[y0,Ynames] = xlsread(file2);

Fpulses = cell(maxtime/1e6*50,1);
Fout = zeros(maxtime/1e6*50,3);
p = 0;
tstart = 0;
tend = 50;
figure(1)
hold on

while tend < tmax
    if rem(tstart/(5*1e5),2) < 1 && rem(tstart,1e4) <50
        tend = tstart+50;
        p = p+1;
        
        PFDon = 300;
        k(1) = PFDon;
    
        if p > 1 
            y0 = Sol.y(:,end);
            
        end
        
        Sol = ode15s(@(t,y)LaiskODENew(t,y,kvals,knames,Ynames),[tstart, tend],y0);
        
%         figure
%         plot(Sol.x, Sol.y)

        Ftmp = Sol.y(end,:); %+ k(fluorescence_reactions(2))*Sol.y(idx2,:);
        ttmp = Sol.x-tstart;
        Fpulses{p} = [Sol.x;Ftmp];
        if ttmp(1) == 0
            ttmp(1) = ttmp(2)/10;
        end
        
%         figure(1)
%         hold on
%         plot(Sol.x,Ftmp)

        logtprime = log(movmean(ttmp,2));
        logtprime = logtprime(2:end);
        dy = diff(Ftmp);
        dx = diff(log(ttmp));
        dydx = dy./dx;
         
        logtprime = logtprime(1:end);
        dydx = dydx(1:end);
         
        idx = (find(ttmp<48,1,'last')-2);
        idxprime = find(logtprime<log(48),1,'last');
 
        pks = findpeaks(dydx(1:idx));
        dvals = sort(pks,'descend');
        maxpeaks = dvals(1:2);
         
        peakidcs = zeros(1:2);
        tpeaks = zeros(1:2);
        peakidcs(1) = find(dydx==maxpeaks(1));
        peakidcs(2) = find(dydx==maxpeaks(2));
         
        dvals2 = sort(peakidcs,'descend');
        furtherdownpeakidc = dvals2(1);
         
        NewTvalues = logtprime(1:furtherdownpeakidc);
        midtimevalue = ((NewTvalues(1)+NewTvalues(end))/(2));
        [~,newtidx] = min(abs(logtprime-midtimevalue));
         
        dydx1=dydx(1:newtidx);
        dydx2=dydx(newtidx:end);
         
        maxpeaksdydx1 = findpeaks(dydx1);
        dvals3 = sort(maxpeaksdydx1,'descend');
        maxpeakleft = dvals3(1:1);
        peakidxleft = find(dydx==maxpeakleft);
         
        maxpeaksdydx2=findpeaks(dydx2);
        dvals4 = sort(maxpeaksdydx2,'descend');
        maxpeakright = dvals4(1:1);
        peakidxright = find(dydx==maxpeakright);
         
        logtmid = ((logtprime(peakidxleft)+logtprime(peakidxright))/(2));
        [~,logtmidx1] = min(abs(logtprime-logtmid));
        tmidx = min(abs(Sol.x - exp(logtprime(logtmidx1))));
        figure(1)
        hold on
%         scatter(tmidx, Ftmp(logtmidx1));
         
        Fo = (Ftmp(logtmidx1));
         
        Fout(p,1) = Fo;
        Fm = Ftmp(find(ttmp<48,1,'last'));
        Fout(p,2) = Fm; 
        Fout(p,3) = (Fm-Fo)/Fm;
        tstart = tend;
 
    elseif rem(tstart/(5*1e5),2) < 1 && rem(tstart,1e4) >= 50
         
        PFD = 0; 
        k(1) = PFD; 
        
        tend = tstart + 1e4 - 50;
        y0 = Sol.y(:,end);
        Sol = ode15s(@(t,y) PSIIODES(t,y,kvals,knames,Ynames),[tstart, tend],y0);
        Ftmp = Sol.y(end,:); %k(fluorescence_reactions(1))*Sol.y(idx1,:);+ k(fluorescence_reactions(2))*Sol.y(idx2,:);
        ttmp = Sol.x-tstart;
        tstart = tend;
%         figure(1)
%         hold on
%         plot(ttmp,Ftmp)


     
    
    elseif rem(tstart/(5*1e5),2) >= 1
        
        PFD = 0;
        k(1) = PFD; 
        
        tend = tstart + 1e4;
        y0 = Sol.y(:,end);
        Sol = ode15s(@(t,y) PSIIODES(t,y,kvals,knames,Ynames),[tstart, tend],y0);
        Ftmp = Sol.y(end,:);%+ k(fluorescence_reactions(2))*Sol.y(idx2,:);
        ttmp = Sol.x-tstart;
        tstart = tend;
%         figure(1)
%         hold on
%         plot(ttmp,Ftmp)

        
    
    end
end
foo = 1;
%     
%     
% for i = 0:1e4:maxtime-1e4
%     
%     to = i;
%     tm = i+50;
%     istart = find(t>=to,1);
%     iend = find(t<=tm,1,'last');
%     if rem(to/(5*1e5),2) < 1 && rem(to,1e4) <50
%         p = p+1; %index of cell array so first pulse will have p = 1, second p = 2
% 
%         ttmp = t(istart:iend)-t(istart);
%         Ftmp = F(istart:iend); %replicate what we have for ttmp and ftmp for everypulse generate a new F0 and Fm
%         Fpulses{p} = [ttmp;Ftmp];
%         if ttmp(1) == 0
%             ttmp(1) = ttmp(2)/10;
%         end
% 
%         logtprime = log(movmean(ttmp,2));
%         logtprime = logtprime(2:end);
%         
%         dy = diff(Ftmp);
%         dx = diff(log(ttmp));
%         dydx = dy./dx;
%         
%         logtprime = logtprime(1:end);
%         dydx = dydx(1:end);
%         
%         idx = (find(ttmp<48,1,'last')-2);
%         idxprime = find(logtprime<log(48),1,'last');
%         %     if isempty(dydx(1:idx))
%         %         pause
%         %     end
%         pks = findpeaks(dydx(1:idx));
%         dvals = sort(pks,'descend');
%         maxpeaks = dvals(1:2);
%         
%         peakidcs = zeros(1:2);
%         tpeaks = zeros(1:2);
%         peakidcs(1) = find(dydx==maxpeaks(1));
%         peakidcs(2) = find(dydx==maxpeaks(2));
%         
%         dvals2 = sort(peakidcs,'descend');
%         furtherdownpeakidc = dvals2(1);
%         
%         NewTvalues = logtprime(1:furtherdownpeakidc);
%         midtimevalue = ((NewTvalues(1)+NewTvalues(end))/(2));
%         [~,newtidx] = min(abs(logtprime-midtimevalue));
%         
%         dydx1=dydx(1:newtidx);
%         dydx2=dydx(newtidx:end);
%         
%         maxpeaksdydx1 = findpeaks(dydx1);
%         dvals3 = sort(maxpeaksdydx1,'descend');
%         maxpeakleft = dvals3(1:1);
%         peakidxleft = find(dydx==maxpeakleft);
%         
%         maxpeaksdydx2=findpeaks(dydx2);
%         dvals4 = sort(maxpeaksdydx2,'descend');
%         maxpeakright = dvals4(1:1);
%         peakidxright = find(dydx==maxpeakright);
%         
%         logtmid = ((logtprime(peakidxleft)+logtprime(peakidxright))/(2));
%         [~,logtmidx1] = min(abs(logtprime-logtmid));
%         scatter(logtprime(logtmidx1), F(logtmidx1));
%         
%         Fo = (Ftmp(logtmidx1));
%         
%         Fout(p,1) = Fo;
%         Fm = Ftmp(find(ttmp<48,1,'last'));
%         Fout(p,2) = Fm; %
%         Fout(p,3) = (Fm-Fo)/Fm;
%         
%         % (log(t(realTvalueidx)),F(realTvalueidx))
%         
%         % figure(1)
%         % hold on
%         % plot(ttmp+(p-1)*100, Ftmp)
%         
%         % scatter([ttmp(io(2)) ttmp(find(ttmp>45,1))]+(p-1)*100 , [Ftmp(io(2)) Ftmp(find(ttmp>45,1))])
%     end
% end
 
 
end 