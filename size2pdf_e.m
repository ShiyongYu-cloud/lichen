function epdfs = size2pdf_e(Site,Xmin,Xmax,dx)
%% function for estimating the empirical pdf of lichen diameters 
% INPUT
% Site: structure containing lichen diameters measured at different sites
% Xmin: lower bound of size domain
% Xmax: upper bound of size domain
% dx: size step
%%
X = Xmin:dx:Xmax;
M = length(X);
N = length(Site);
epdfs = zeros(M,N);
%% get results
for i = 1:N
    prob = epdf(Site(i).size,Xmin,Xmax,dx);
    prob(prob<0.0001) = 0; 
    epdfs(:,i) = prob;
end
return;
%%
function prob = epdf(x,Xmin,Xmax,dx)
%% function for estimating empirical pdf 
%INPUT
%x: vector specifying lichen sizes measured at a site 
%Xmin: lower bound of the size domain
%Xmax: upper bound of the size domain
%dx: spatial resolution 
%OUTPUT
%prob: vector specifying estimated PDF at spatial points
%%
% generating spatial points at which pdf will be estimated
X = Xmin:dx:Xmax;
% estimating the CDF by definition
M = length(X);
N = length(x);
prob = zeros(M,1);
F = zeros(M,1);
h = 0.0000001;        % Step for evaluating thenumerical differentiation
for i = 1:M
    p = 0;              % True Probability
    q = 0;              % False Probability
    for j = 1:N
        if x(j) <= X(i)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(i) = p/(p + q);   % Calulating Probability
end
% Estimating PDF by differentiating the CDF
for k = 1:M
    fxph = interp1(X,F,X(k) + h,'spline');  % Interpolating value of F(x+h)
    fxmh = interp1(X,F,X(k) - h,'spline');  % Interpolating value of F(x-h)
    prob(k) = (fxph - fxmh)/(2*h); 
    if prob(k) < 0
       prob(k) = 0;
    end        
end                                         
% normalizing to 1
prob = smooth(prob);  
prob = prob./sum(prob);
return;