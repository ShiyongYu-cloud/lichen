function par_hat = parest(Xmin,Xmax,Tmin,Tmax,dx,dt,F,S,par_0,par_min,par_max)
%% function for estimating model parameters
%INPUT
% F: matrix containing the PDFs of the observational data
% par_0: initial value of parameters
% par_min: minimum value of parameters
% par_max: maximum value of parameters
%OUTPUT 
% par_hat = optimal estimate of parameters
%% 
myopt = optimset('fminsearch');
myopt = optimset(myopt,'MaxFunEvals',50000,'MaxIter',5000,'TolFun',1.e-3,'TolX',1.e-3);
par_hat = fminsearchbnd(@(par) costfun(par,Xmin,Xmax,Tmin,Tmax,dx,dt,F,S),par_0,par_min,par_max,myopt);
end