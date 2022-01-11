function Div = costfun(par,Xmin,Xmax,Tmin,Tmax,dx,dt,F,Site)
%% function for calculating the total Jensen-Shannon divergence of two PDFs
%INPUT
% par: vector containing model parameters [a,b,sigma]
% Xmin: lower bound of the size domain
% Xmax: upper bound of the size domain
% Tmin: lower bound of the time domain
% Tmax: upper bound of the time domain
% dx: size step
% dt: time step
% F: matrix containing the PDFs of lichen diameters measured at different sites
% Site: structure containing the lichen diameters measured at different sites
%OUTPUT
% Div: sum of the Jensen-Shannon divergences
%% setting up computing grids
T = Tmin:dt:Tmax; % temporal grid points
%% solving the PDE of lichen growth 
P = growthpde(par,Xmin,Xmax,Tmin,Tmax,dx,dt);
%% calculating the total Jensen-Shannon divergence
M = length(Site); % number of study sites
t = zeros(1,M);
for i = 1:M 
    t(i) = 2000 - Site(i).age;
end
% find the nearest grid points of t in T
t_T = round(t/dt)*dt; % round ages to the nearest dt years
dist = ones(1,M);
ind = ismember(T,t_T);
P_hat = P(:,ind);
for i = 1:M 
    dist(i) = D_JS(P_hat(:,i),F(:,i));
end    
Div = sum(dist);
end