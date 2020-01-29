function [FvFm, grads] = calc_FvFm(x0,... Set of parameters. This only includes the independent variables as described by the third column in Y and Constants files
    n_trains, n_flashes, flash_duration, flash_interval, train_interval, ... Experimental parameters
    Fluorescence_k_idcs, Fluorescence_y_inds,... indeces used to calculate fluorescence
    kidcs, PSIidcs, ... all indices needed in to calculate FvFm and prepare the variables
    tablek, tabley,... information on the k and y variables
    kconst, rate_inds, S, species, knames, species_idcs) % model specific variables



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



k(kidcs.PFD) = .24;
k(kidcs.kf) = 1;
PS1T = k(kidcs.a2)*k(kidcs.Chl)/k(kidcs.PSU2);
PS2T = (1-k(kidcs.a2))*(k(kidcs.Chl)/k(kidcs.PSU1));
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
grads = zeros(n_trains*n_flashes,1);
mult1 = kidcs.mult1;
mult2 = kidcs.mult2;
n1idx = kidcs.n1idx;
% dark adapt the system
k(mult1) = 0;
k(mult2) = 0;
k(n1idx) = 0;
dark_adaptation_time = 180; %3-5 minutes typically
t_lims = [0,dark_adaptation_time];
Sol =  ode15s(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t_lims,yinitial);
yinitial = Sol.y(:,end); %initialize the y vector for the next iteration
counter = 1;
for train = 1:n_trains
    %     fprintf('train %i \n', train)
    for flash = 1:n_flashes
        k(mult1) = mult1Val;
        k(mult2) = mult2Val;
        k(n1idx) = n1;
        nTimepoints = flash_duration*1e6;
        t = linspace(0, flash_duration, nTimepoints);
        Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
        while any(any(Sol<-1e-5)) || any(any(isnan(Sol)))
            nTimepoints = nTimepoints*5;
            t = linspace(0, flash_duration, nTimepoints);
            Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
            
        end
        Sol = Sol';
        F = LaiskFluorescence(Fluorescence_y_inds, Fluorescence_k_idcs, k, Sol) ;
        FvFm(counter) = (max(F) - F(1))/max(F);
        grads(counter) = mean(abs(diff(F(end-5:end))./diff(t(end-5:end))));
        counter = counter+1;
        yinitial = Sol(:,end); %initialize the y vector for the next iteration
        %Shift to dark time between flashes
        k(mult1) = 0;
        k(mult2) = 0;
        k(n1idx) = 0;
        t_lims = [0,flash_interval];
        Sol = ode15s(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t_lims,yinitial);
        if any(any(Sol.y<-1e-5)) || any(any(isnan(Sol.y)))
            nTimepoints = flash_interval*1e3;
            t = linspace(0, flash_interval, nTimepoints);
            Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
            while any(any(Sol<-1e-5)) || any(any(isnan(Sol)))
                nTimepoints = nTimepoints*5;
                t = linspace(0, flash_interval, nTimepoints);
                Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
            end
            Sol = Sol';
            yinitial = Sol(:,end);
        else
            yinitial = Sol.y(:,end);
        end
        
    end
    k(mult1) = 0;
    k(mult2) = 0;
    k(n1idx) = 0;
    t = linspace(0, train_interval, train_interval*1e5);
    Sol = ode15s(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
    if any(any(Sol.y<-1e-5)) || any(any(isnan(Sol.y)))
        nTimepoints = train_interval*1e3;
        t = linspace(0, train_interval, nTimepoints);
        Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
        while any(any(Sol<-1e-5)) || any(any(isnan(Sol)))
            nTimepoints = nTimepoints*5;
            t = linspace(0, train_interval, nTimepoints);
            Sol = ode2(@(t,y) PS2ODES(t,y,k(kconst),k,rate_inds,S,knames,species),t,yinitial);
        end
        Sol = Sol';
        yinitial = Sol(:,end);
    else
        yinitial = Sol.y(:,end);
    end
    %     yinitial = Sol.y(:,end); %initialize the y vector for the next iteration
end
end

% figure; plot(1:length(FvFm), FvFm, 'o-')
% foo = 1;