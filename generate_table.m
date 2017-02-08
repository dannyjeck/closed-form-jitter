function table = generate_table(D,Nmax)
% table = generate_table(D,Nmax) 
% Generates a table of probabilities P(C|Nx,Ny,D) for a fixed value of D and
% an optional maximum value of Nx and Ny. C can range from 0 to D.
%
% INPUTS: 
% D    - the jitter window length in number of bins. Usually each bin is 1 ms.
% Nmax - restricts the size of the output table to be [D+1 x Nmax+1 x Nmax+1].
%        Detault value is D.
%
% OUTPUTS:
% table - table of probabilities. table(c+1,nx+1,ny+1) = P(c|nx,ny,D)

if nargin<2
    Nmax = D;
end

%Generate table of P(C|Nx,Ny,D).
table = zeros(D+1,Nmax+1,Nmax+1); %use D+1 because number of coincidences can be 0
for ny = 0:Nmax
    for nx = 0:Nmax
        cmax = min(nx,ny); %maximum possible number of coincidences in this window
        cmin = max(0,nx-(D-ny));
        for c = cmin:cmax
            table(c+1,nx+1,ny+1) = nchoosek(D-ny,nx-c)*nchoosek(ny,c);
        end
    end
end
table = table./repmat(sum(table),[D+1 1 1]);