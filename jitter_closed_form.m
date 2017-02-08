function [excess_synch, p] = jitter_closed_form(X,Y,D,taumax,table,table_mean)
% [excess_synch, p] = jitter_exact(X,Y,D,taumax,table,table_mean) computes the
% jitter-corrected correlogram and performs a significance test on from two signals 
% at a given jitter time scale D. taumax,table, and table_mean are optional
% parameters.
%
% INPUTS:
% X,Y     - two binary spike trains (1 at the time bins where a spike
%           occurs, 0 otherwise).
% D       - time scale of rate covariations to be tested
% taumax  - maximum lag to be tested. outputs will be 
%           created for each value in the range of (-taumax:taumax). 
%           A positive lag means X was shifted to the left.
% table   - A pre-computed table of probabilites P(C_int|Nx,Ny,D) assuming the
%           null hypothesis (that spike times can be jittered within a window) 
%           is true. If not provided, the table will be computed within
%           this function.
% table_mean- A pre-computed table of expected values E[P(C_int|Nx,Ny,D)] assuming 
%           the null hypothesis (that spike times can be jittered within a window) 
%           is true. If not provided, table_mean will be computed within
%           this function. 
%
% OUTPUTS:
% excess_synch - the number of synchronous spikes beyond the mean 
%       of the null hypothesis distribution. 
% p -   a vector of p-values testing whether the measured correlation is significantly
%       different from the null hypothesis distribution at each value of -taumax:taumax.
%       If the significance test is not requested, processing is much faster.
%
% Author: Daniel Jeck
% Copywrite: March, 2014

%Variable Glossary:
% Nwin   - Number of jitter intervals (sometimes called windows)
% Nx, Ny - spike counts in X or Y for the various jitter intervals 
% Nmax   - Maximum value of Nx or Ny
% count  - matrix fo the number of times (Nx,Ny) occurs in the signal
% C_true - measured correlation at different lags
% Cmax   - maximum possible correlation value
% speclen- maximum correlation value in the final probability distribution
% pdf    - distribution of the number of coincidences under the jitter null hypothesis

%truncate to the nearest window length and make column vectors
if D<2
    error('Jitter window D must be at least 2 bins');
end

Nwin = floor(length(X)/D);
if Nwin<1
    error('Signal length must be greater than one jitter interval');
end
X = X(1:D*Nwin); X = double(X(:));
Y = Y(1:D*Nwin); Y = double(Y(:));

if nargin<4 %taumax not given, choose the length of the smallest signal as taumax
    taumax = min(length(X),length(Y));
end
if taumax>length(X)
    error('taumax too big for input signal');
end
% Compute actual correlation function
intflag = all(abs([X; Y]-round([X; Y])<eps));
C_true = xcorr(double(X),double(Y),taumax);
if intflag
    C_true = round(C_true);
end
tau_list = -taumax:taumax;

%Compute Ny for all intervals
Ny = sum(reshape(Y,D,Nwin));

Nmax = max(max(conv(X,ones(D,1))),max(Ny)); %convolution of X gives all possible values of Nx for all values of tau
Cmax = min(sum(X),sum(Y)); %max possible number of coincidences in the whole signal

%Generate table of P(C_int|Nx,Ny,D)
%Note about the tables: They include values for Nx=0, Ny=0, and C_int = 0.
%Therefore, they will be indexed with +1's where necessary.
if nargin<5
    table = generate_table(D,Nmax);
else
    table = table(:,1:Nmax+1,1:Nmax+1);
end

if nargout<2
    if nargin<6
        %compute expected value for each value of (Nx,Ny)
        w = repmat((0:D)',[1 Nmax+1 Nmax+1]); 
        table_mean = squeeze(sum(table.*w)); 
    else
        table_mean = table_mean(1:Nmax+1,1:Nmax+1);
    end
end


p = ones(length(tau_list),1);
excess_synch = ones(length(tau_list),1);
for t = 1:length(tau_list)
    
    %compute Nx for shifted version of X
    if tau_list(t)<0
        X_shift = [zeros(-tau_list(t),1); X(1:end+tau_list(t))];
    else 
        X_shift = [X(1+tau_list(t):end); zeros(tau_list(t),1)];
    end
    Nx = sum(reshape(X_shift,D,Nwin));
    
    %Count incidices of the values in the table;
    count = zeros(Nmax+1,Nmax+1);
    for j = 1:Nwin
        count(Nx(j)+1,Ny(j)+1) = count(Nx(j)+1,Ny(j)+1) +1;
    end
    
    
    if nargout<2
        mean_H0 = sum(table_mean(:).*count(:)); %mean number of coincidences under the null hypothesis
        excess_synch(t) = C_true(t)-mean_H0;
    else
        %pdf of the number of coincidences for the full signal (under the null hypothesis) 
        % is the convolution the interval pdfs. This convolution is done in
        % the frequency domain for speed.
        
        % Take padded FFTs of useful rows, apply exponent, and multiply
        speclen = max(D,Cmax);
        table_pad = [table; zeros(speclen-D,Nmax+1,Nmax+1)];
        table_pad = table_pad(:,:);
        count_list = count(:)>0;
        spec = fft(table_pad(:,count_list));
        spec = spec.^repmat(count(count_list)',speclen+1,1);
        spec = prod(spec,2);

        %take iFFT to get closed form pdf
        pdf = ifft(spec);
        pdf = pdf/sum(pdf); %to correct for numerical errors

        p(t) = sum(pdf(C_true(t)+1:Cmax+1));
        excess_synch(t) = C_true(t)-(0:speclen)*pdf;
    end
end
% 
