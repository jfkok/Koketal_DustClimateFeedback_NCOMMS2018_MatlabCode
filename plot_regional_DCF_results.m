%delete(findall(0,'Type','figure'));

%plotting the regional DCF shown in Fig. 5 in paper
figure(1); clf;
map_plot_script_zero_log_subplot(1000*DCF_CMIP5median_tot_median,10^2.5,1,1,2,1,[2,2,1],'RdBu',''); 
set(gca,'FontSize',20); %title('Dust-climate feedback based on CMIP5 median (mW m^{-2} K^{-1})','FontSize',26)
xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
map_plot_script_zero_log_subplot(1000*DCF_CMIP5lowCI_tot_median,10^2.5,1,1,2,1,[2,2,2],'RdBu','');
set(gca,'FontSize',20); %title('Dust-climate feedback based on CMIP5 lower CI (mW m^{-2} K^{-1})','FontSize',26)
xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
map_plot_script_zero_log_subplot(1000*DCF_CMIP5upCI_tot_median,10^2.5,1,1,2,1,[2,2,3],'RdBu','');
set(gca,'FontSize',20); %title('Dust-climate feedback based on CMIP5 upper CI (mW m^{-2} K^{-1})','FontSize',26)
xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
map_plot_script_zero_log_subplot(1000*DCF_CESMnew_tot_median,10^2.5,1,1,2,1,[2,2,4],'RdBu','');
set(gca,'FontSize',20); %title('Dust-climate feedback based on CESM with PHYS module (mW m^{-2} K^{-1})','FontSize',26)
xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');

%plotting the SW and LW DRE based on each model, Figure S5 in paper
k=1; %CESM
map_plot_script_zero_log(squeeze(DRE_SW_corr_avg(k,:,:)),10,3,2,1,k+4,'RdBu','SW dust DRE based on CESM (W/m^2)'); set(gca,'FontSize',32);
map_plot_script_zero_log(squeeze(DRE_LW_corr_avg(k,:,:)),10,3,2,1,k+no_DRE_models+4,'RdBu','LW dust DRE based on CESM (W/m^2)'); set(gca,'FontSize',32);
k=2; %GISS
map_plot_script_zero_log(squeeze(DRE_SW_corr_avg(k,:,:)),10,3,2,1,k+4,'RdBu','SW dust DRE based on GISS (W/m^2)'); set(gca,'FontSize',32);
map_plot_script_zero_log(squeeze(DRE_LW_corr_avg(k,:,:)),10,3,2,1,k+no_DRE_models+4,'RdBu','LW dust DRE based on GISS (W/m^2)'); set(gca,'FontSize',32);
k=3; %GEOS-Chem
map_plot_script_zero_log(squeeze(DRE_SW_corr_avg(k,:,:)),10,3,2,1,k+4,'RdBu','SW dust DRE based on GEOS-Chem (W/m^2)'); set(gca,'FontSize',32);
map_plot_script_zero_log(squeeze(DRE_LW_corr_avg(k,:,:)),10,3,2,1,k+no_DRE_models+4,'RdBu','LW dust DRE based on GEOS-Chem (W/m^2)'); set(gca,'FontSize',32);
k=4; %WRF-Chem
map_plot_script_zero_log(squeeze(DRE_SW_corr_avg(k,:,:)),10,3,2,1,k+4,'RdBu','SW dust DRE based on WRF-Chem (W/m^2)'); set(gca,'FontSize',32);
map_plot_script_zero_log(squeeze(DRE_LW_corr_avg(k,:,:)),10,3,2,1,k+no_DRE_models+4,'RdBu','LW dust DRE based on WRF-Chem (W/m^2)'); set(gca,'FontSize',32);

%plotting the regional DRE: total, SW, and LW; Figure S6 in paper
map_plot_script_zero_log(DRE_regional_SW_median,10,3,2,1,5,'RdBu','SW dust DRE (W/m^2)'); set(gca,'FontSize',36);
map_plot_script_zero_log(DRE_regional_LW_median,10,3,2,1,6,'RdBu','LW dust DRE (W/m^2)'); set(gca,'FontSize',36);
map_plot_script_zero_log(DRE_regional_tot_median,10,3,2,1,7,'RdBu','Total dust DRE (W/m^2)'); set(gca,'FontSize',36);

%plotting the regional DCF enhancement: total, SW, and LW;
%Figure S7 in the paper
map_plot_script_zero_log(DCF_enhancement_SW_median,100,2,2,2,8,'PRGn','Regional enhancement of SW dust climate feedback'); set(gca,'FontSize',32);
map_plot_script_zero_log(DCF_enhancement_LW_median,100,2,2,2,9,'PRGn','Regional enhancement of LW dust climate feedback'); set(gca,'FontSize',32);
map_plot_script_zero_log(DCF_enhancement_tot_median,10^1.5,1,1,2,10,'PRGn','Regional enhancement of total dust climate feedback'); set(gca,'FontSize',32);

%plotting the SW and LW regional DCF using CESMnew
map_plot_script_zero_log(1000*DCF_CESMnew_SW_median,10^2.5,1,1,2,11,'RdBu','SW DCF in mW m^{-2} K^{-1} based on CESM');
map_plot_script_zero_log(1000*DCF_CESMnew_LW_median,10^2.5,1,1,2,12,'RdBu','LW DCF in mW m^{-2} K^{-1} based on CESM');

%plotting the SW and LW regional DCF using CMIP5
map_plot_script_zero_log(1000*DCF_CMIP5median_SW_median,10^2.5,1,1,2,13,'RdBu','SW DCF in mW m^{-2} K^{-1} based on CMIP5');
map_plot_script_zero_log(1000*DCF_CMIP5median_LW_median,10^2.5,1,1,2,14,'RdBu','LW DCF in mW m^{-2} K^{-1} based on CMIP5');

    

%outdated code:
% %the regional total dust DRE:
% figure(1); clf; map_plot_script_vdf (lon,lat,DRE_regional_tot_median,-10,10,'Total DRE in W/m');
% 
% 
% 
% axis_lim = 5;
% DRE_regional_LW_median_plot=DRE_regional_LW_median;
% DRE_regional_LW_median_plot(DRE_regional_LW_median_plot<-axis_lim)=-axis_lim;
% DRE_regional_LW_median_plot(DRE_regional_LW_median_plot>axis_lim)=axis_lim;
% figure(3); clf; map_plot_script_vdf_relative([lon-180;180],lat,([DRE_regional_LW_median_plot;DRE_regional_LW_median_plot(73,:)]),-axis_lim,axis_lim,'LW DRE in W/m');
% axis_lim = 5;
% DRE_regional_SW_median_plot=DRE_regional_SW_median;
% DRE_regional_SW_median_plot(DRE_regional_SW_median<-axis_lim)=-0.999*axis_lim;
% DRE_regional_SW_median_plot(DRE_regional_SW_median>axis_lim)=axis_lim;
% figure(4); clf; map_plot_script_vdf_relative([lon-180;180],lat,([DRE_regional_SW_median_plot;DRE_regional_SW_median_plot(73,:)]),-axis_lim,axis_lim,'SW DRE in W/m');
% 
% %plotting the regional DCF enhancement of the total DRE:
% figure(3); clf; map_plot_script_vdf (lon,lat,DCF_enhancement_tot_median,-40,40,'Regional enhancement of global DCF)');
% plot_dust_AOD (xlon_plot([73:144,1:72]),lat,DCF_enhancement_tot_median([73:144,1:72],:),mean(mean(DCF_enhancement_tot_median)),5,'Regional DCF enhancement',3);
% plot_dust_AOD (xlon_plot([73:144,1:72]),lat,abs(DCF_enhancement_tot_median([73:144,1:72],:)),mean(mean(abs(DCF_enhancement_tot_median))),5,'Regional DCF enhancement',4);
% 
% %plotting the results:
% DRE_LW_corr_avg_plot = squeeze(DRE_LW_corr_avg);
% plot_dust_AOD (xlon_plot([73:144,1:72]),lat,DRE_LW_corr_avg_plot([73:144,1:72],:),mean(mean(DRE_LW_corr_avg_plot)),5,'LW DRE (W/m2)',1);
% 
% DRE_SW_corr_avg_plot = squeeze(DRE_SW_corr_avg);
% figure(2); clf; map_plot_script_vdf(lon,lat,DRE_SW_corr_avg_plot,-5,5,'SW DRE (W/m2)');
% 
% DRE_tot_corr_avg_plot = DRE_LW_corr_avg_plot+DRE_SW_corr_avg_plot;
% figure(3); clf; map_plot_script_vdf (lon,lat,DRE_tot_corr_avg_plot,-5,5,'Tot DRE (W/m2)');
% plot_dust_AOD (xlon_plot([73:144,1:72]),lat,DRE_tot_corr_avg_plot([73:144,1:72],:),mean(mean(DRE_tot_corr_avg_plot)),1,'Tot DRE (W/m2)',3);
% 
% DRE_tot_corr_norm_avg_plot = squeeze(DRE_tot_corr_norm_avg);
% figure(4); clf; map_plot_script_vdf (lon,lat,DRE_tot_corr_norm_avg_plot,-200,200,'Tot DRE, normalized by total global DRE (W/m2)');
