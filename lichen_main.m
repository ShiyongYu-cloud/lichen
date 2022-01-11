%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inverse Model of Lichen Growth Curves                                 %%
%% Copyright: Shiyong Yu                                                 %% 
%% E-mail: shiyong.yu@gmail.com                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc;
%% Entering the data for sites of known age (calibration data)
SizeUnit = 'mm'; %Unit of lichen size
StudyArea = 'Huashan'; %Study area
AgeScale = 'b2k'; % Year scale of the ages to be reported (b2k or CE)
delta = 5; % nearest ages to be rounded 
Site(1).age = 1984;
Site(2).age = 1976;
Site(3).age = 1961;
Site(4).age = 1844;
Site(5).age = 1835;
Site(6).age = 1820;
Site(7).age = 1787;
Site(8).age = 1776;
Site(9).age = 1762;
Site(10).age = 1696;
Site(11).age = 1670;
Site(12).age = 1593;
Site(1).size = [20,23,25,20,24,30,26,26,17,19,23];
Site(2).size = [26,29,32,29,31,32,27,29,22,35];
Site(3).size = [33,36,41,26,44,38,31,33,35,34];
Site(4).size = [70,63,61,56,66,64,58,63];
Site(5).size = [64,66,67,68,65,65,76,61,71];
Site(6).size = [68,74,72,79,70,68,65,71,73,69];
Site(7).size = [82,80,78,75,80,85,82,80,83];
Site(8).size = [79,85,92,85,83,75,77,86,96,85,68];
Site(9).size = [84,87,81,75,86,81,91,94];
Site(10).size = [101,96,99,102,95,104,105,101,100,93,102,93,107,102,101,98,89,94,110];
Site(11).size = [109,102,104,106,98,99,102,105,113];
Site(12).size = [116,118,123,122,125,124,127,134];
%% Entering the data for sites to be datad
sample(1).ID = 'Grotto No. 24';
sample(2).ID = 'Grotto No. 2';
sample(1).size = [92,101,82,90,85,87,91,84,86];
sample(2).size = [41,55,47,48,51,56,44];
%% Specifying the computing grids for the growth curve
dx = 1;         %spatial resolution of 1 cm)
dt = 1;         %temporal resolution of 1 year
%%
Xmin = 0;       
Tmin = 0;
M = length(Site);
MaxSize = zeros(1,M);
MeanSize = zeros(1,M);
for i = 1:M
   MaxSize(i) = max(Site(i).size);
   MeanSize(i) = mean(Site(i).size);
end
Xmax = ceil(max(MaxSize)/100)*100;
t = zeros(1,M);
% convert ages to the b2k scale
for i = 1:M 
    t(i) = 2000 - Site(i).age;
end
t = round(t/dt)*dt; % round to the nearest dt years
Tmax = ceil(t(end)/100)*100;
%% Building up the PDF of the lichen size datasets
F = size2pdf_gev(Site,Xmin,Xmax,dx);
%% Estimating model parameters of lichen growth
y = diff(MeanSize)./diff(t);
y = y';
x = [1./sqrt(t(1:end-1)') ones(M-1,1)];
beta0 = x\y; 
sigma0 = sqrt(var(y-x*beta0));
par_0 = [beta0',sigma0];
par_min = [0,0,0];
par_max = [10,1,1];
par_hat = parest(Xmin,Xmax,Tmin,Tmax,dx,dt,F,Site,par_0,par_min,par_max);
%% Calculating and plotting the growth curve
%dx = 1; 
%dt = 1;
X = Xmin:dx:Xmax; % size grid points
X = X';
T = Tmin:dt:Tmax;
P = growthpde(par_hat,Xmin,Xmax,Tmin,Tmax,dx,dt);
figure(1);
plot_curve(X,T,P,Site,Xmax,Tmax,SizeUnit,AgeScale);
%% Estimating the empirical PDF of ages of the undated sites
F_sample = size2pdf_e(sample,Xmin,Xmax,dx);
P(isnan(P)) = 0;
N = length(T);
m = length(sample);
prob = zeros(N,m);
for i = 1:N
    for j = 1:m
        prob(i,j) = 1 - D_JS(P(:,i),F_sample(:,j));
    end
end 
prob(prob < 0.001) = 0;
prob = prob./(ones(N,1)*sum(prob));
%% Calculating the HPD intervals of ages and plotting the results
AGES = pdf2hpd(T,prob,delta,sample,AgeScale,StudyArea);
figure(2);
plot_ages(T,prob,AGES,AgeScale);
%% Cleaning up
clear beta0 delta dt dx F F_samples i j m M MaxSize MeanSize N par_0 par_min par_max ...
    sigma0 t Tmax Tmin x Xmax Xmin y F_sample