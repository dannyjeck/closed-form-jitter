function [excess_synch, p] = jitter_monte_carlo(X,Y,D,Ntrials,taumax)
% [excess_synch, p] = jitter_monte_carlo(X,Y,D,Ntrials,taumax)
% computes the jitter-corrected correllogram (JCCG) and performs a 
% significance test on from two signals at a given jitter time scale D. 
%
% INPUTS:
% X,Y     - two binary spike trains (1 at the time bins where a spike
%           occurs). A positive lag means X has been shifted left
% D       - time scale of rate covariations to be tested
% Ntrials - number of surrogate spike trains to be generated
% taumax  - maximum lag to be tested. outputs will be 
%           created for each value in the range of (-taumax:taumax). 
%           A positive value of taumax means X was shifted to the left.
%
% OUTPUTS:
% excess_synch - the number of synchronous spikes beyond the mean 
%       of the null hypothesis distribution. 
% p -   a p-value testing whether the measured synchrony is significantly
%       different from the null hypothesis distribution. 


if nargin<5 %if taumax not given, use the length of the shortest signal
    taumax = min(length(X),length(Y));
end

%truncate to the nearest window length and make column vectors
Nwin = floor(length(X)/D);
X = X(1:D*Nwin); X = double(X(:));
Y = Y(1:D*Nwin); Y = double(Y(:));

% Compute actual correlation function
intflag = all(abs([X; Y]-round([X; Y])<eps));
C_true = xcorr(double(X),double(Y),taumax);
if intflag
    C_true = round(C_true);
end

%generate surrogate spike trains of X by jittering in each interval
X_MC = zeros(length(X),Ntrials);
for j = 1:Nwin 
    idx = (j-1)*D+1:j*D;
    X_win = X(idx);
    
    tmp = arrayfun(@(x)randperm(D),(1:Ntrials)','UniformOutput',0); %generate random permutations
    X_MC(idx,:) = X_win(cell2mat(tmp)');
    
end

% Process differently if you only need a portion of the xcorr (because xcorr.m
% computes the whole cross-correlation and then outputs a fraction of it)
if nargin<4
    C_trials = cell2mat(arrayfun(@(n)xcorr(Y,X_MC(:,n),taumax),1:Ntrials,'UniformOutput',0));
else
    X_MC = [zeros(taumax,Ntrials); double(X_MC); zeros(taumax,Ntrials)];
    C_trials = conv2(X_MC,flipud(double(Y)),'valid');
end
% 

%compute p values for all lags simultaneously
p = (sum(repmat(C_true,1,Ntrials)<=C_trials,2)+1)/(Ntrials+1);

%compute JCCG estimate
excess_synch = C_true-mean(C_trials,2);
