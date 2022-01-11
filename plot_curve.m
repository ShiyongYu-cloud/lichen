function [] = plot_curve(X,T,P,Site,Xmax,Tmax,SizeUnit,AgeScale)
%% Function for plotting the growth curve
%INPUT
%X: size grids
%T: time grids
%Site: structure containing the measurements of lichen sizes
%Xmax: upper bound of spatial domain
%Tmax: upper bound of temporal domain
%SizeUnit: unit of lichen size (mm or cm)
%AgeScale: scale of age to be reported (b2k or CE)
%%
P(P<0.001) = NaN;
if strcmpi(AgeScale,'b2k') == 1
    pcolor(T,X,P);
    colormap(hsv);
    shading interp;
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca, 'TickDir', 'out');
    ylabel(strcat('Largest thallus diameter', ' (', SizeUnit, ')'));
    xlabel(['Surface exposure age ', '(years b2k)']);
    ylim([0 Xmax]);
    xlim([0 Tmax]);
    h = colorbar;
    caxis([0 ceil(max(max(P(:,10:end)))/0.1)*0.1]);
    ylabel(h,'Probability');
    grid on;
    %% Plotting the measurements
    M = length(Site);
    for i = 1:M
        n = length(Site(i).size);
        hold on;
        plot(2000-Site(i).age*ones(1,n),Site(i).size,'o','MarkerEdgeColor','k','MarkerFaceColor','w'); 
    end
elseif strcmpi(AgeScale,'CE') == 1
    T = 2000 - T;
    pcolor(T,X,P);
    colormap(hsv);
    shading interp;
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca, 'TickDir', 'out');
    set(gca,'Xdir','reverse');
    ylabel(strcat('Largest thallus diameter', ' (', SizeUnit, ')'));
    xlabel(['Surface exposure age ', '(CE)']);
    ylim([0 Xmax]);
    xlim([2000-Tmax 2000]);
    h = colorbar;
    caxis([0 ceil(max(max(P(:,10:end)))/0.1)*0.1]);
    ylabel(h,'Probability');
    grid on;
    %% Plotting the measurements
    M = length(Site);
    for i = 1:M
        n = length(Site(i).size);
        hold on;
        plot(Site(i).age*ones(1,n),Site(i).size,'o','MarkerEdgeColor','k','MarkerFaceColor','w'); 
    end    
end
end