% Speed_test.m
% By: Daniel Jeck
% Copywrite October 2013
%
% This script runs a simulation of the type that was used to generate
% Figure 1 in "Closed Form Jitter Analysis of Neuronal Spike Trains" 
% by Daniel Jeck and Ernst Niebur (http://arxiv.org/abs/1502.07907).
% It compares the time it takes to perform a Monte Carlo Jitter
% analysis with the Closed Form Jitter Method in computing P values

T_list = [1 1 1 1 1 1];% 1:10:100; %list of singal lengths (In seconds) to be tested
f_list = [2 10]; %[5 10 20 40 50 100 200 400 500]; %list of firing frequencies (in Hz) to be tested
repeat = 2; %number of signals tested per condition

times_closed_form = zeros(length(f_list),length(T_list),repeat); %matrix of computation times for closed form algorithm
times_jitter_mc= zeros(length(f_list),length(T_list),repeat); %matrix of computation times for Monte Carlo algorithm
for tidx = 1:length(T_list)
    for fidx = 1:length(f_list);
        for rep = 1:repeat
            %% Generate signals
            dt = .001; %time bin size (seconds)
            T = T_list(tidx);
            t = 0:dt:T-dt;
            f = f_list(fidx);
            p_fire = dt*f;
            X1 = rand(size(t))<p_fire;
            X2 = X1;
            
            %% Jitter (closed form and Monte Carlo) parameters
            taumax = 100; %maximum correlation lag of interest (in time bins)
            D = 20; %jitter interval width (in time bins)
            Ntrials = 1000; %Number of Monte Carlo Trials

            %% Measure computation time
            tic
                [excess2, p2] = jitter_closed_form(X1,X2,D,taumax);
            times_closed_form(fidx,tidx,rep) = toc;

            tic
                [excess, p] = jitter_monte_carlo(X1,X2,D,Ntrials,taumax);

            times_jitter_mc(fidx,tidx,rep) = toc;
        end
        display([fidx tidx]);
    end
end
%% Average over the repeats
times_jitter_mc_mean = mean(times_jitter_mc,3);
times_closed_form_mean = mean(times_closed_form,3);

%Correct for the fact that 1000 MC runs were done when 20,000 were needed
times_jitter_mc_mean_bf = 20*times_jitter_mc_mean; 

%% Find cases where MC jitter is faster
compare = times_jitter_mc_mean_bf<=times_closed_form_mean;
line = nan(1,length(T_list));
for k = 1:length(T_list)
    a = find(compare(:,k),1,'first');
    if ~isempty(a)
        line(k) = f_list(a);
    end
end

%% Absolute time plots (not in paper)
figure(1);
Imax = max(max(max(times_jitter_mc_mean_bf)),max(max(times_closed_form_mean)));
subplot(121);
imagesc(T_list,f_list,times_closed_form_mean,[0 Imax])
hold all;plot(T_list,line,'k');hold off;
xlabel('Spike Train Length (s)');
ylabel('Spike frequency (Hz)');
title('Closed Form Jitter Method');
axis square

subplot(122);
imagesc(T_list,f_list,times_jitter_mc_mean_bf,[0 Imax]);colorbar
hold all;plot(T_list,line,'k');hold off;
xlabel('Spike Train Length (s)');
% ylabel('Spike frequency (Hz)');
title('Monte Carlo Jitter');
axis square

%% Speedup ratio plot
figure(2);
semilogy(T_list,(times_jitter_mc_mean_bf(:,:)./times_closed_form_mean(:,:))','k')
axis square;
for k = 1:length(f_list)
    f = f_list(k);
    txt = [num2str(f) ' Hz'];
    h_txt(k) = text(T_list(6),(times_jitter_mc_mean_bf(k,6)./times_closed_form_mean(k,6)),txt,'BackgroundColor','w');
    
end
xlabel('Spike Train Length [s]');
ylabel('Performance Gain');
title('P value Computation');

%find minimum and maximum speed up ratios
disp([min((times_jitter_mc_mean_bf(:)./times_closed_form_mean(:))) max((times_jitter_mc_mean_bf(:)./times_closed_form_mean(:)))]);

