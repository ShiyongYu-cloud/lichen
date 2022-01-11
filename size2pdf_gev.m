function pdfs = size2pdf_gev(Site,Xmin,Xmax,dx)
%% function for estimating the GEV distribution of lichen diameters 
% INPUT
% Site: structure containing lichen diameters measured at different sites
% Xmin: lower bound of the lichen diameter
% Xmax: upper bound of the lichen diameter
% dx: size step
%%
X = Xmin:dx:Xmax;
X = X';
M = length(X);
N = length(Site);
pdfs = zeros(M,N);
MU = zeros(N,1);
S = zeros(N,1);
K = zeros(N,1);
%% obtain gev parameters
for i = 1:N
    dia = Site(i).size;
    dia = dia';
    pd = fitdist(dia,'gev');
    MU(i) = pd.mu;
    S(i) = pd.sigma;
    K(i) = pd.k;
end
%% obtain pdf at the grid points
for i = 1:N
    pdfs(:,i) = gevpdf(X,K(i),S(i),MU(i));
end
end