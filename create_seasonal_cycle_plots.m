%plots the seasonal cycle in AOD at each AERONET station

m = ceil(sqrt(total_stations));
n = ceil(total_stations/m);

for p=1:total_stations
    s = seasonal_station(p);
    if (s==1)
    	h=figure(2); clf; 
	end
    clear mod_x_plot; clear mod_y_plot; clear obs_plot; clear seas_cycle_plot;
    obs_plot(1:j)=station_monthly_save(1,s,:); %obs_plot(isnan(obs_plot))=-1;
    obs_plot(obs_plot==0)=NaN;
    x_plot = 1:j;
    figure(2); subplot(m,n,p); hold on; plot(x_plot, obs_plot,'k-','LineWidth',1); %errorbar(obs_plot,obs_error,'ok','LineWidth',1);
    y_max = max(obs_plot);
    y_min = min(obs_plot);
    seas_cycle_plot(1:12)=avg_monthly_AOT(s,:); 
    seas_cycle_plot=repmat(seas_cycle_plot,1,year_end-year_start+1); %extending the (12 point) seasonal cycle to the full data range
    plot(x_plot,seas_cycle_plot,'r','LineWidth',0.5);
    seas_cycle_corr_plot=obs_plot-seas_cycle_plot;
    plot(x_plot,seas_cycle_corr_plot,'g','LineWidth',0.5);
    y_max = 1.02*max(obs_plot);
    y_min = 1.02*min(seas_cycle_corr_plot);
    start_j=min(find(~isnan(obs_plot)));
    axis([start_j,j,1.25*y_min,1.25*y_max]);
    set(gca,'FontSize',9);
    text(1.1,1.18*y_max,'r = ','FontSize',8); %bold font to indicate significance
    if (station_lat(s)>0) %in NH
        string = strcat(station_names(s,:),', ',num2str(station_lat(s)),'N,  ',num2str(station_lon(s)),'W');
    else %in SH
        string = strcat(station_names(s,:),', ',num2str(-station_lat(s)),'S,  ',num2str(station_lon(s)),'W');
    end
    title(string,'FontSize',9,'FontWeight','bold');
    xlabel('Month','FontSize',9);
    ylabel('AERONET AOD','FontSize',9);
end %for
