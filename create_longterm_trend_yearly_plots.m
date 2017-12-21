%this plots the long-term records of AOD at the dusty stations, as well as
%their anomalies after subtracting the seasonal cycle

m = ceil(sqrt(total_longterm_stations));
n = ceil(total_longterm_stations/m);
no_years = year_end-year_start+1;

for p=1:total_longterm_stations
    s = has_longterm_trend(p);
    %plotting the observed and modeled AOD
    if (p==1)
    	h = figure(1); clf; figure(2); clf;
    end
    %calculate annually-averaged AOD for obs and model runs from seasonal cycle
    obs_avg_AOD(p) = sum([31,28,31,30,31,30,31,31,30,31,30,31].*AOT_seas_cycle(has_longterm_trend(q),:))/365;
    model_avg_AOD(1,p) = sum([31,28,31,30,31,30,31,31,30,31,30,31]'.*squeeze(model_TotAOT_seas_cycle(1,has_longterm_trend(p),:)))/365;
    model_avg_AOD(2,p) = sum([31,28,31,30,31,30,31,31,30,31,30,31]'.*squeeze(model_TotAOT_seas_cycle(2,has_longterm_trend(p),:)))/365;
    %clearing variables for re-use
    clear mod_x_plot; clear mod_y_plot; clear obs_plot; 
    clear mod_y2_plot; clear obs2_plot;
    %for plot 1: anomolies
    x_plot = year_start:year_end;
    obs_plot=meas_annual_mean_AOT_seas_corr(s,:); %obs_plot(isnan(obs_plot))=-1;
    figure(1); h=subplot(m,n,p); hold on; plot(x_plot, obs_plot,'ko-','LineWidth',1); %errorbar(obs_plot,obs_error,'ok','LineWidth',1);
    figure(1); hold on; plot(x_plot,AER_intercept(s)+AER_slope(s)*(x_plot-x_plot(1)+1),'k--');
    y_max1 = max(obs_plot);
    y_min1 = min(obs_plot);
    %for plot 1: absolute values
    obs_plot2=obs_avg_AOD(p)+meas_annual_mean_AOT_seas_corr(s,:); %obs_plot(isnan(obs_plot))=-1;
    figure(2); h=subplot(m,n,p); hold on; plot(x_plot, obs_plot2,'ko-','LineWidth',1); %errorbar(obs_plot,obs_error,'ok','LineWidth',1);
    %figure(2); hold on; plot(x_plot,AER_intercept(s)+AER_slope(s)*(x_plot-x_plot(1)+1),'k--');
    y_max2 = max(obs_plot);
    y_min2 = min(obs_plot);
    for r=1:no_model_runs
        mod_y_plot(1:size(x_plot,2))=model_annual_mean_AOT_seas_corr(r,s,:); 
        mod_y2_plot(1:size(x_plot,2))=model_avg_AOD(r,p)+model_annual_mean_AOT_seas_corr(r,s,:); 
        if (r==1)
            figure(1); plot(x_plot,mod_y_plot,'ro-','LineWidth',0.5);
            figure(1); plot(x_plot,mod_intercept(r,s)+mod_slope(r,s)*(x_plot-x_plot(1)+1),'r--');
            figure(2); plot(x_plot,mod_y2_plot,'ro-','LineWidth',0.5);
        elseif (r==2)
            figure(1); plot(x_plot,mod_y_plot,'mo-','LineWidth',0.5);
            figure(1); plot(x_plot,mod_intercept(r,s)+mod_slope(r,s)*(x_plot-x_plot(1)+1),'m--');
            figure(2); plot(x_plot,mod_y2_plot,'mo-','LineWidth',0.5);
        elseif (r==3)
            figure(1); plot(x_plot,mod_y_plot,'go-','LineWidth',0.5);
            figure(1); plot(x_plot,mod_intercept(r,s)+mod_slope(r,s)*(x_plot-x_plot(1)+1),'g--');
            figure(2); plot(x_plot,mod_y2_plot,'go-','LineWidth',0.5);
        elseif (r==4)
            figure(1); plot(x_plot,mod_y_plot,'bo-','LineWidth',0.5);
            figure(1); plot(x_plot,mod_intercept(r,s)+mod_slope(r,s)*(x_plot-x_plot(1)+1),'b--');
            figure(2); plot(x_plot,mod_y2_plot,'bo-','LineWidth',0.5);
        elseif (r==5)
            figure(1); plot(x_plot,mod_y_plot,'co-','LineWidth',0.5);
        end
        if (max(mod_y_plot)>y_max1)
            y_max1 = double(max(mod_y_plot));
        end
        if (min(mod_y_plot)<y_min1)
            y_min1 = double(min(mod_y_plot));
        end        
        if (max(mod_y2_plot)>y_max2)
            y_max2 = double(max(mod_y2_plot));
        end
        if (min(mod_y2_plot)<y_min2)
            y_min2 = double(min(mod_y2_plot));
        end        
    end
    start_j=min(find(~isnan(obs_plot)))+year_start-1;
    figure(1); axis([year_start,year_end,-0.15,0.15]);
    set(gca,'FontSize',9);
    string = strcat('Fraction of variance explained by Angstrom exponent = ',num2str(AOD_angstrom_corr(s)^2),','); %string for writing out correlation between measured AOD and Angstrom exponent
    text(year_end-0.62*(year_end-start_j),1.05*y_min1,string,'FontSize',8); %writing out correlation between measured AOD and Angstrom exponent
    text(start_j+0.02*(year_end-start_j),1.12*y_max1,'r = ','FontSize',8); %writing out correlation between modeled and measured AOD
    for r=1:no_model_runs
        switch r
            case 1
                string=strcat(num2str(trend_corrcoef(r,s),'%2.2f'),',');
                figure(1); text(start_j+0.12*(year_end-start_j),1.12*y_max1,string,'FontSize',8,'Color','r','FontWeight','bold'); 
            case 2
                string=strcat(num2str(trend_corrcoef(r,s),'%2.2f'),',');
                figure(1); text(start_j+0.22*(year_end-start_j),1.12*y_max1,string,'FontSize',8,'Color','m','FontWeight','bold'); 
            case 3
                string=strcat(num2str(trend_corrcoef(r,s),'%2.2f'),',');
                figure(1); text(start_j+0.32*(year_end-start_j),1.12*y_max1,string,'FontSize',8,'Color','g','FontWeight','bold');
            case 4
                string=num2str(trend_corrcoef(r,s),'%2.2f');
                figure(1); text(start_j+0.42*(year_end-start_j),1.12*y_max1,string,'FontSize',8,'Color','b','FontWeight','bold');
        end
    end
    if (station_lat(s)>0 && station_lon(s)<180) %in NH and east
        string = strcat(station_names(s,:),', ',num2str(station_lat(s)),'N,  ',num2str(station_lon(s)),'E');
    elseif (station_lat(s)>0 && station_lon(s)>=180) %in NH and west
        string = strcat(station_names(s,:),', ',num2str(station_lat(s)),'N,  ',num2str(360-station_lon(s)),'W');
    elseif (station_lat(s)<0 && station_lon(s)>=180) %in SH and east
        string = strcat(station_names(s,:),', ',num2str(-station_lat(s)),'S,  ',num2str(station_lon(s)),'E');
    else %in SH and west
        string = strcat(station_names(s,:),', ',num2str(-station_lat(s)),'S,  ',num2str(360-station_lon(s)),'W');
    end
    title(string,'FontSize',9,'FontWeight','bold');
    xlabel('Month','FontSize',9);
    ylabel('Aeos','FontSize',9);
    
    %plotting the observed angstrom exponent
    clear mod_x_plot; clear mod_y_plot; clear obs_plot;
    obs_plot=meas_annual_mean_angstrom(s,:); %obs_plot(isnan(obs_plot))=-1;
    obs_error=zeros(1,no_years); obs_error(obs_error==0)=NaN;
    x_plot = year_start:year_end;
    h=figure(1); h=subplot(m,n,p); hold on; plot(x_plot, obs_plot,'ko-','LineWidth',1); %errorbar(obs_plot,obs_error,'ok','LineWidth',1);
    figure(1); hold on; plot(x_plot,angstrom_intercept(s)+angstrom_slope(s)*(x_plot-x_plot(1)+1),'k--');
end %for

text(start_j+1.12*(year_end-start_j),0.75*y_max1,'Average r = ','FontSize',16,'FontWeight','bold'); %bold font to indicate significance
for r=1:no_model_runs
    switch r
    	case 1
        	string=strcat(num2str(mean(trend_corrcoef(r,has_longterm_trend)),'%2.2f'),',');
            text(start_j+1.62*(year_end-start_j),0.75*y_max1,string,'FontSize',13,'Color','r','FontWeight','bold'); %bold font to indicate significance
        case 2
        	string=strcat(num2str(mean(trend_corrcoef(r,has_longterm_trend)),'%2.2f'),',');
            text(start_j+1.82*(year_end-start_j),0.75*y_max1,string,'FontSize',13,'Color','m','FontWeight','bold'); %bold font to indicate significance
        case 3
            string=strcat(num2str(mean(trend_corrcoef(r,has_longterm_trend)),'%2.2f'),',');
            text(start_j+2.02*(year_end-start_j),0.75*y_max1,string,'FontSize',13,'Color','g','FontWeight','bold'); %bold font to indicate significance
        case 4
        	string=num2str(mean(trend_corrcoef(r,has_longterm_trend)),'%2.2f');
            text(start_j+2.22*(year_end-start_j),0.75*y_max1,string,'FontSize',13,'Color','b','FontWeight','bold'); %bold font to indicate significance
	end
end
1;