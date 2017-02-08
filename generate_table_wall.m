function table = generate_table_wall(D,w,Nmax)
%Generates a table of probabilities P(C|Nx,Ny) for a fixed value of D and
%an optional maximum value of Nx and Ny. C can range from 0 to D.
%
% INPUTS: 
% D - the jitter window length in number of bins. Usually each bin is 1 ms.
% w - the odds ratio of a coincidence when there are an even number of
%     possibilities. A Wallenius noncentral hypergeometric distribution is
%     used.
% Nmax - restricts the size of the output table to be [D+1 x Nmax+1 x Nmax+1].
%        Detault value is D.
%
% OUTPUTS:
% table - table of probabilities. table(c+1,nx+1,ny+1) = P(c|nx,ny)

if nargin<3
    Nmax = D;
end

%Generate table of P_i(C|Nx,Ny)
table = zeros(D+1,Nmax+1,Nmax+1); %use D+1 because number of coincidences can be 0
table(1,1,:) = 1;
for ny = 0:Nmax
    for nx = 1:Nmax
        cmax = min(nx,ny); %maximum possible number of coincidences in this window
        cmin = max(0,nx-(D-ny));
        
        %% Note: may be room for improvement to vectorize this
        for c = cmin:cmax
            if c ==0
                table(c+1,nx+1,ny+1) = table(c+1,nx,ny+1).* (D-ny-nx+1+c)./((D-ny-nx+1+c)+w*(ny-c));
            else
                table(c+1,nx+1,ny+1) = table(c+1,nx,ny+1).* (D-ny-nx+1+c)./((D-ny-nx+1+c)+w*(ny-c))...
                                     + table(c  ,nx,ny+1).* w.*(ny-c+1)./(D-ny-nx+c+w.*(ny-c+1));
            end
        end
    end
end

%normalizing may be redundant
table = table./repmat(sum(table),[D+1 1 1]);
