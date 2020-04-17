function [ts,ys, Fs, FvFm, O2, O2_light, O2_dark] = calc_Species_concs_ode15s1(x0,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
                    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
                    Fluorescence_k_idcs, Fluorescence_y_inds,... indeces used to calculate fluorescence
                    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
                    tablek, tabley,... information on the k and y variables
                    kconst, rate_inds, S, species, knames, species_idcs, Rknames, analysis_name,yidcs,ATPpar,kf1indcs, kf2indcs) % model specific variables
                

indepy = find(tabley.independent);
% yinitial = zeros(length(species),1);
yinitial = tabley.base_val;
yr = x0(1:length(indepy));
yinitial(indepy) = yr;
yinitial = yinitial(species_idcs);

indepk = find(tablek.independent);
kr = x0(length(indepy)+1:end);
k = tablek.base_val;
k(indepk) = kr;


% k(kidcs.PFD) = .015;
k(kidcs.kf) = 1;
PS1T = (1-k(kidcs.a2))*k(kidcs.Chl)/k(kidcs.PSU1); 
PS2T = k(kidcs.a2)*(k(kidcs.Chl)/k(kidcs.PSU2));
n1 = k(kidcs.PFD)*k(kidcs.Labs)*(1-k(kidcs.a2))/PS1T;
n2 = k(kidcs.PFD)*k(kidcs.Labs)*k(kidcs.a2)/PS2T;                 
yinitial(PSIidcs) = PS1T/PS2T;
     
mult1Val = n2*k(kidcs.kp)/(1+k(kidcs.kp)+k(kidcs.kn)+k(kidcs.kr));
mult2Val = n2*k(kidcs.kp)/(1+k(kidcs.kp)+k(kidcs.kn));
Div1Val = k(kidcs.kpc)/k(kidcs.kEpc);
Div2Val = k(kidcs.kcytf)/k(kidcs.kEcytf);
Div3Val = k(kidcs.kfx)/k(kidcs.kEfx);
Div4Val = k(kidcs.kb6f)/k(kidcs.kEb6f);

k(kidcs.mult1) = mult1Val;
k(kidcs.mult2) = mult2Val;
k(kidcs.Div1) = Div1Val;
k(kidcs.Div2) = Div2Val;
k(kidcs.Div3) = Div3Val; 
k(kidcs.Div4) = Div4Val;
k(kidcs.n1idx) = n1;
FvFm = zeros(n_trains*n_flashes,1);

mult1 = kidcs.mult1;
mult2 = kidcs.mult2;
n1idx = kidcs.n1idx;
ts = {};
ys = {};
Fs = {};
rs = {};
O2 =  zeros(n_trains*n_flashes,1);
O2_light =  zeros(n_trains*n_flashes,1);
O2_dark =  zeros(n_trains*n_flashes,1);
O2ind = find(strcmp(species, 'O2'));
% dark adapt the system
k(mult1) = 0;
k(mult2) = 0;    
k(n1idx) = 0;
dark_adaptation_time = 1; %3-5 minutes typically
t_lims = [0,dark_adaptation_time];
yinitial(yidcs.deltapsiindex) = 0;
yinitial(yidcs.ATPaseoindex) = .63;
yinitial(yidcs.ATPaserindex) = 0; 
yinitial(yidcs.pH_lumenindex) = 7.8;
yinitial(yidcs.pH_stromaindex) = 7.8;
yinitial(yidcs.fRindex) = 0;
species{end+1} = 'deltapsi';
species{end+1} = 'ATPaseO';
species{end+1} = 'ATPaseR';
species{end+1} = 'pHLumen';
species{end+1} = 'pHStroma';
species{end+1} = 'fr';

% Sol =  ode15s(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t_lims,yinitial);
% ts{end+1} = -dark_adaptation_time+Sol.x;
% ys{end+1} = Sol.y;
% Fs{end+1} = [];
ts{end+1} = 0;%-dark_adaptation_time+Sol.x;
ys{end+1} = 0;%Sol.y;
Fs{end+1} = [];
% yinitial = Sol.y(:,end); %initialize the y vector for the next iteration 
counter = 1;
for train = 1:n_trains
    fprintf('train %i \n', train)
    if train == 9
        foo = 1;
    end
    for flash = 1:n_flashes
        fprintf('flash %i \n', flash)
        k(mult1) = mult1Val;
        k(mult2) = mult2Val;    
        k(n1idx) = n1;
        nTimepoints = flash_duration*1e6;
        t = [0, flash_duration];
        Sol = ode15s(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);
        while any(any(Sol.y<-1e-5)) || any(any(isnan(Sol.y)))
            nTimepoints = nTimepoints*5;
            t = linspace(0, flash_duration, nTimepoints);
            Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);    
            
        end
        ts{end+1} = ts{end}(end) + Sol.x;

        ys{end+1} = Sol.y;

% 
%         t = linspace(0, flash_duration, nTimepoints);
%         Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);
%         ts{end+1} = ts{end}(end) + t;
% 
%         Sol = Sol';
%         ys{end+1} = Sol';

        if length(ts{end}) ~= size(ys{end},2)
            foo = 1;
        end
        F = LaiskFluorescence(Fluorescence_y_inds, Fluorescence_k_idcs, k, Sol.y) ;
%         F = LaiskFluorescence(Fluorescence_y_inds, Fluorescence_k_idcs, k, Sol) ;
        F1 = F;
        t1 = Sol.x;
%         t1 = t;
        FvFm(counter) = (max(F) - F(1))/max(F);
        flash_O2 = Sol.y(O2ind,:) -Sol.y(O2ind,1) ;
%         flash_O2 = Sol(O2ind,:) -Sol(O2ind,1) ;
        O2(counter) = trapz(Sol.x, flash_O2);
%         O2(counter) = trapz(t, flash_O2);
        O2_light(counter) = trapz(Sol.x, flash_O2);
%         O2_light(counter) = trapz(t, flash_O2);
        yinitial = Sol.y(:,end); %initialize the y vector for the next iteration 
        Fs{end+1} = F;
        if any(isnan(yinitial))
            foo = 1;
        end

        %Shift to dark time between flashes
        k(mult1) = 0;
        k(mult2) = 0;      
        k(n1idx) = 0;
        
        t_lims = [0,flash_interval];
        Sol = ode15s(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t_lims,yinitial);
%         yinitial = Sol.y(:,end); %initialize the y vector for the next iteration 
        if any(any(Sol.y<-1e-5)) || any(any(isnan(Sol.y)))
            nTimepoints = flash_interval*1e3;
            t = linspace(0, flash_interval, nTimepoints); 
            Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);    
            while any(any(Sol<-1e-5)) || any(any(isnan(Sol)))
                nTimepoints = nTimepoints*5;
                t = linspace(0, flash_interval, nTimepoints);
                Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);    
            end
            Sol = Sol';
            ts{end+1} = ts{end}(end) + t;
            %         ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+t;
            ys{end+1} = Sol;
            Fs{end+1} = [];
            yinitial = Sol(:,end);
        else
            ts{end+1} = ts{end}(end) + Sol.x;
            %         ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+Sol.x;
            ys{end+1} = Sol.y;
            Fs{end+1} = [];
            yinitial = Sol.y(:,end);
        end



%         nTimepoints = flash_interval*1e5;
%         t = linspace(0, flash_interval, nTimepoints); 
%         Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);
%         ts{end+1} = ts{end}(end) + t;
%         Fs{end+1} = [];
%         Sol = Sol';
%         yinitial = Sol(:,end);
%         O2_dark_prod = Sol(O2ind,:) - Sol(O2ind,1);
%         O2_dark(counter) = trapz(t, O2_dark_prod);
%         O2(counter) = O2(counter)+O2_dark(counter);
        
%         
%         O2_dark_prod = Sol.y(O2ind,:) - Sol.y(O2ind,1);
%         O2_dark(counter) = trapz(Sol.x, O2_dark_prod);
%         O2(counter) = O2(counter)+O2_dark(counter);
        counter = counter+1;
        if length(ts{end}) ~= size(ys{end},2)
            foo = 1;
        end
    end
    k(mult1) = 0;
    k(mult2) = 0;       
    k(n1idx) = 0;
    t = linspace(0, train_interval, train_interval*1e5);
    Sol = ode15s(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);
%     ts{end+1} = Sol.x;
        if any(any(Sol.y<-1e-5)) || any(any(isnan(Sol.y)))
                nTimepoints = train_interval*1e3;
                t = linspace(0, train_interval, nTimepoints);
                Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);    
                while any(any(Sol<-1e-5)) || any(any(isnan(Sol)))
                    nTimepoints = nTimepoints*5;
                    t = linspace(0, train_interval, nTimepoints);
                    Sol = ode2(@(t,y) PS2ODES1(t,y,k(kconst),k,rate_inds,S,Rknames,species,yidcs,ATPpar,kf1indcs, kf2indcs,kidcs),t,yinitial);
                end
            Sol = Sol';
            ts{end+1} = ts{end}(end) + t;
            %         ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+t;
            ys{end+1} = Sol;
            Fs{end+1} = [];
            yinitial = Sol(:,end);
        else
            ts{end+1} = ts{end}(end) + Sol.x;
            %         ts{end+1} = (train-1)*(train_interval + n_flashes*(flash_duration+flash_interval)) + n_flashes*(flash_duration+flash_interval)+Sol.x;
            ys{end+1} = Sol.y;
            Fs{end+1} = [];
            yinitial = Sol.y(:,end);
        end
%     if length(ts{end}) ~= size(ys{end},2)
%             foo = 1;
%     end
%     ys{end+1} = Sol.y; %calculate the species evolution during the dark between flashes
%     Fs{end+1} = []; % save an empty vector for fluorescence (to allign the values of times and fluorescence)
     %initialize the y vector for the next iteration    
end
disp('Computations Complete. Plotting ...')

% 
% figure; plot(1:length(O2_dark), O2_dark, '.-');
% title('O2 produced in the dark')
% 
% figure; plot(1:length(O2_dark), O2_light, '.-');
% title('O2 produced in the light')
% 
% 
% figure;plot(1:length(FvFm), FvFm, '.-')
% ylabel('Flash FqFm')
% yyaxis right;
% plot(1:length(O2), O2, '.-')
% ylabel('Flash O_2 yield')
% xlabel('STF #')
% title('Combined O_2 and FvFm oscillations')
% 

Ty = removevars(tabley, {'lb', 'ub', 'independent'});
Tk = removevars(tablek, {'lb', 'ub', 'independent'});
titles = {};

figure;
semilogx(t1,F1)
ylabel('Fluorescence')
xlabel('time')
titles{end+1} = 'Fluorescence';
savefig('Fluorescence.fig')

figure;
plot(1:length(FvFm), FvFm, '.-')
ylabel('Flash FqFm')
xlabel('STF #')
title('FvFm oscillations')
titles{end+1} = 'FvFm_oscillations';
savefig('FvFm_oscillations.fig')

figure;
plot(1:length(O2), O2, '.-')
ylabel('Flash O_2 yield')
xlabel('STF #')
title('O_2 oscillations')
titles{end+1} = 'O2_oscillations';
savefig('O2_oscillations.fig')

figure;
plot(1:length(FvFm), FvFm, '.-')
ylabel('Flash FqFm')
yyaxis right
plot(1:length(O2), O2, '.-')
ylabel('Flash O_2 yield')
xlabel('STF #')
title('Combined FvFm and O2 oscillations')
titles{end+1} = 'FvFm_O2_oscillations';
savefig('Combined FvFm_O2_oscillations.fig')

plot_S_states(species, ys, ts);
titles{end+1} = 'S_states';
plot_pq_redox_state(species, ys, ts);
titles{end+1} = 'PQ';
plot_cyt_redox_state(species, ys, ts);
titles{end+1} = 'Cyt';
plot_pc_redox_state(species, ys, ts);
titles{end+1} = 'PC';
plot_P700_redox_state(species, ys, ts);
titles{end+1} = 'P700';
plot_fd_redox_state(species, ys, ts);
titles{end+1} = 'Fd';
plot_H_species(species, ys, ts);
titles{end+1} = 'H_species';
plot_NAD_redox_state(species, ys, ts);
titles{end+1} = 'NADP';
plot_atp_species(species, ys, ts);
titles{end+1} = 'ADP_phospho';
titles{end+1} = 'ADP_ATP';
plot_H_species_ATPSYN(species, ys, ts);
titles{end+1} = 'H_species_ATP_SYN';
plot_H2O_species(species, ys, ts)
titles{end+1} = 'water_species';
plot_buffer_species(species, ys, ts)
titles{end+1} = 'buffer_species';
plot_OH_species(species, ys, ts)
titles{end+1} = 'OH_species';
plot_Fl(Fs, ts);
fm = reshape(FvFm,n_flashes,[]);
a = mean(fm,1);
figure; plot(1:length(a), a,'.-')
ylabel('Average Fq/Fm per train')
xlabel('train number')
titles{end+1} = 'Avg_FvFm';
titles{end+1} = 'dsadas';


figure;
plot(t1,F1);
ylabel('Fluorescence')
xlabel('time')
titles{end+1} = 'First_Flash_Fluorescence';


h = get(0,'children');
timestamp = datestr(now,30);

if ~exist('results', 'dir')
    mkdir('results')
end
if ~exist(['results/' analysis_name], 'dir')
    mkdir('results/', analysis_name);
end

mkdir(['results/', analysis_name, '/', timestamp]);
dirname = ['results/', analysis_name, '/', timestamp];

for i=1:length(h)
    saveas(h(i), [dirname, '/', titles{i}], 'png');
%     saveas(h(i), [dirname, '/', titles{i}], 'fig');
end

writetable(Ty, [dirname, '/ys.csv'])
writetable(Tk, [dirname, '/ks.csv'])


end

% figure; plot(1:length(FvFm), FvFm, 'o-')
% foo = 1;