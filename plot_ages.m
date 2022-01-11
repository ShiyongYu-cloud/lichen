function [] = plot_ages(T,prob,AGES,AgeScale)
%% function for plotting the calendar ages
%INPUT
%T: time grids
%prob: Jensen-Shannon similarity
%AGES: structure containing the estimated ages of the samples
%AgeScale: scale of age to be reported (BP or CE)
%%
if strcmpi(AgeScale,'CE') == 1
    T = 2000 - T;
end
cal_age = T';
N = size(prob,2);
maxprob = max(max(prob))*N/(N-1);
for i = 1:N
    prob(:,i) = prob(:,i) + maxprob*(i-1);  % blow up and shift
    ind1 = (prob(:,i) > maxprob*(i-1));
    early = max(cal_age(ind1));
    late = min(cal_age(ind1));
    ind2 = (cal_age <= early) & (cal_age >= late);
    calage = cal_age(ind2);
    CALAGE = [calage(1); calage; calage(end); calage(1)];
    PROB = [maxprob*(i-1); prob(ind2,i); maxprob*(i-1); maxprob*(i-1)];
    fill(CALAGE,PROB,[0.301 0.745 0.933],'EdgeColor','none');  % plot pdfs
    hold on
    p95_4 = AGES(i).P95_4_Credible_intervals;    % plot the 95.4% pdf
    K = size(p95_4,1);
    for j = 1:K
        id = (cal_age >= p95_4(j,1)) & (cal_age <= p95_4(j,2));
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob(id,i);
        PDF = [maxprob*(i-1); pdf; maxprob*(i-1); maxprob*(i-1)];
        fill(AGE,PDF,[0.929 0.694 0.125],'EdgeColor','none');
    end
    hold on
    p68_2 = AGES(i).P68_2_Credible_intervals;    % plot the 68.2% pdf 
    L = size(p68_2,1);
    for k = 1:L
        id = (cal_age >= p68_2(k,1)) & (cal_age <= p68_2(k,2));
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf = prob(id,i);
        PDF = [maxprob*(i-1); pdf; maxprob*(i-1); maxprob*(i-1)];
        fill(AGE,PDF,[0.635 0.078 0.184],'EdgeColor','none'); 
    end
    line(calage, prob(ind2,i),'Color','k','LineWidth',0.1);
    line(T,prob(:,i),'Color','k','LineWidth',0.1);
end
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickDir', 'out');
ylim([0 N*maxprob]);
set(gca,'yticklabel',[]);
set(gca,'ytick',[]);
if strcmpi(AgeScale,'b2k') == 1 
    xlim([min(cal_age) max(cal_age)]);
    set(gca,'XDir','reverse');
    xlabel(['Surface exposure age ', '(years b2k)']);
elseif strcmpi(AgeScale,'CE') == 1
    xlim([round(min(T)/100)*100 round(max(T)/100)*100]);
    xlabel('Surface exposure age (CE)');  
end
% plot sampleID 
for i = 1:N
    if strcmpi(AgeScale,'b2k') == 1 
       text(max(cal_age)-25,maxprob*0.8+maxprob*(i-1),AGES(i).Sample_ID); 
    elseif strcmpi(AgeScale,'CE') == 1 
       text(min(cal_age)+25,maxprob*0.8+maxprob*(i-1),AGES(i).Sample_ID); 
    end
end
end