function AGES = pdf2hpd(T,prob,delta,sample,AgeScale,StudyArea) 
%% function for calculating the HPD intervals
%INPUT: 
%T: time grids
%prob: Jensen-Shannon similarity at the time grids
%delta: nearest years to be rounded
%sample: structure containing the size of samples to be dated
%AgeScale: scale of year to be reported (b2k or CE)
%OUTPUT
%AGES: structure containing 68.2% and 95.4% CI and median probability age  
%%
M = length(sample);  %total number of samples
AGES = struct('Sample_ID',[],'P68_2_Credible_intervals',[],'P95_4_Credible_intervals',[],'Median_prob_age',[],'Year_scale',[]);
for i = 1:M
    AGES(i).Sample_ID = strcat(StudyArea,{' '},sample(i).ID);
    [AGES(i).P68_2_Credible_intervals,AGES(i).P95_4_Credible_intervals,AGES(i).Median_prob_age] = hpd(T,prob(:,i),delta);
end
if strcmpi(AgeScale,'b2k') == 1
    for i = 1:M
        AGES(i).Year_scale = strcat('Years before',{' '},num2str(2000), {' '},'CE');
        AGES(i).P68_2_Credible_intervals = fliplr(AGES(i).P68_2_Credible_intervals);
        AGES(i).P95_4_Credible_intervals = fliplr(AGES(i).P95_4_Credible_intervals);
        AGES(i).Median_prob_age = fliplr(AGES(i).Median_prob_age);
    end
elseif strcmpi(AgeScale,'CE') == 1
    for i = 1:M
        AGES(i).Year_scale = AgeScale;
        AGES(i).P68_2_Credible_intervals = 2000 - AGES(i).P68_2_Credible_intervals;
        AGES(i).P95_4_Credible_intervals = 2000 - AGES(i).P95_4_Credible_intervals;
        AGES(i).Median_prob_age = 2000 - AGES(i).Median_prob_age;
    end
end
return;
%%
function [p68_2,p95_4,median_age] = hpd(T,prob,delta)
%% function for estimating HPD intervals  
%INPUT
%T: time grids
%prob: JS similarity at the time grids   
%delta: nearest age to be rounded
%OUTPUT
%p68_2: 68.2% calendar age probability
%p95_4: 95.4% calendar age probability
%median_age: calendar age at median probability
%%
% Code is adopted from MatCal 3.0 (2020-08-13)
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%% Calculating HPD intervals
T = T(:);
prob = prob(:);
calprob = [T prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
%1 sigma credible interval
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);
hpd68_2 = sortrows(hpd68_2,1);
ind1 = find(diff(hpd68_2(:,1)) > delta);
if isempty(ind1) == 1
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
	%p68_2(1,3) = sum(hpd68_2(1:end,2));
    %p68_2(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind1)
		z = z + 1;
		indy1(z,1) = ind1(i);
		z = z + 1;
		indy1(z,1) = ind1(i)+1;
	end
	indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
	z=0;
	for i = 2:2:length(indy1)
		z = z+1;
		p68_2(z,1) = hpd68_2(indy1(i),1);
		p68_2(z,2) = hpd68_2(indy1(i-1),1);
		%p68_2(z,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));
	end
	p68_2 = flipud(p68_2);
    %p68_2(:,3) = p68_2(:,3)/sum(p68_2(:,3)); %fraction of each enclosed area in the pdf
end
%2 sigma credible interval
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);
hpd95_4 = sortrows(hpd95_4,1);
ind2 = find(diff(hpd95_4(:,1)) > delta);
if isempty(ind2) == 1
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
	%p95_4(1,3) = sum(hpd95_4(1:end,2));
    %p95_4(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind2)
		z = z + 1;
		indy2(z,1) = ind2(i);
		z = z + 1;
		indy2(z,1) = ind2(i)+1;
	end
	indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
	z=0;
	for i = 2:2:length(indy2)
		z = z+1;
		p95_4(z,1) = hpd95_4(indy2(i),1);
		p95_4(z,2) = hpd95_4(indy2(i-1),1);
		%p95_4(z,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));
	end
	p95_4 = flipud(p95_4);
    %p95_4(:,3) = p95_4(:,3)/sum(p95_4(:,3)); %fraction of enclosed area in the pdf
end
%% calculate age at median probability
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
median_age = round(calprob(median_ind(1),1));
return;