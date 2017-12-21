clear all;

%this function reads in the different CLM runs, extracts the global dust
%flux, and calibrates it using the calibration to AERONET data in Kok et
%al. (2014)

start_year = 2000;
end_year = 2100;
running_avg_years = 25;
model_dir = '..\Data\CESM\CLM output\'; %location of model output

pi = 3.14159;
r_earth = 6.371*10^6; %earth radius in meters
year_in_seconds = 60*60*24*365;

case_old_param = 'qian48_mb95_fecan'; %with old dust emission parameterization and Fecan et al. (1999) soil moisture effects
case_new_param = 'qian48_kok_fecan'; %with new dust emission parameterization and Fecan et al. (1999) soil moisture effects
path_old_param = strcat(model_dir,case_old_param,'\');
path_new_param = strcat(model_dir,case_new_param,'\');

%Global dust emission rates previously calculated from calibration of CESM/CAM4/CLM4-CN runs in Kok et al. (2014) from ERA-I runs (1995-2011):
calib_global_dust_flux_old_param = 1e3*4.0257; %global dust emission rate in ERA-I run in Kok et al. (2014), for original parameterization with new soil moisture, in Tg/year
calib_global_dust_flux_new_param = 1e3*3.1753; %global dust emission rate in ERA-I run in Kok et al. (2014), for new parameterization with new soil moisture, in Tg/year

%reading in the area for the CLM output grid boxes
ncid = netcdf.open(strcat(model_dir,'surfdata_1.9x2.5_simyr1850_c091108.nc'),'nowrite');
varid_area = netcdf.inqVarID(ncid,'AREA'); %in km2
area = 1e6*netcdf.getVar(ncid,varid_area); %the 1e6 converts from km2 to m2
netcdf.close(ncid);

%reading in the source function for 2.5 degree lon x 1.9 degree lat
ncid = netcdf.open(strcat(model_dir,'dst_1.9x2.5_c090203.nc'),'nowrite');
varid_source_fct = netcdf.inqVarID(ncid,'mbl_bsn_fct_geo');
source_fct = netcdf.getVar(ncid,varid_source_fct);
netcdf.close(ncid);

%reading in forcing files lat and lon, and calculating the area of each gridbox
ncid = netcdf.open(strcat(model_dir,'clmforc.Qian.c2006.T62.TPQW.2000-01.nc'),'nowrite');
varid_lon = netcdf.inqVarID(ncid,'LONGXY');
lon = netcdf.getVar(ncid,varid_lon);
varid_lat = netcdf.inqVarID(ncid,'LATIXY');
lat = netcdf.getVar(ncid,varid_lat);
netcdf.close(ncid);

delta_lat = 180/(size(lat,2)); %grid box extent in latitudinal degrees
delta_lon = 360/(size(lon,1)); %grid box extent in longitudinal degrees
delta_lat_length = 2*pi*r_earth*delta_lat/360; %meridional grid box extent in meters 
delta_lon_length = 2*pi*r_earth*delta_lon/360; %zonal grid box extent in meters

for i=1:size(lon,1)
    for j=1:size(lat,2)
        temperature_area(i,j) = cos(pi*lat(i,j)/180)*delta_lat_length*delta_lon_length;
	end
end

for p=1:end_year-start_year+running_avg_years
    year(p) = start_year + p - running_avg_years; year(p)
        
    %calculating global dust flux for case old param
    filename_old_param = strcat(path_old_param,'dust_ea_',case_old_param,'_',num2str(year(p)),'.nc');
    ncid = netcdf.open(filename_old_param,'nowrite');
    varid_dust = netcdf.inqVarID(ncid,'DSTFLXT');
    dust_flux = netcdf.getVar(ncid,varid_dust);
    dust_flux(dust_flux > 1e35) = 0; %getting rid of the fill values and replacing them with zeroes
    global_dust_emis_rate_old_param(p) = sum(sum(dust_flux.*area.*source_fct*year_in_seconds))/10^9; %global dust emission rate in Tg/year; old parameterization uses source function

    %calculating global dust flux for case new param
    filename_new_param = strcat(path_new_param,'dust_ea_',case_new_param,'_',num2str(year(p)),'.nc');
    ncid = netcdf.open(filename_new_param,'nowrite');
    varid_dust = netcdf.inqVarID(ncid,'DSTFLXT');
    dust_flux = netcdf.getVar(ncid,varid_dust);
    dust_flux(dust_flux > 1e35) = 0; %getting rid of the fill values and replacing them with zeroes
    netcdf.close(ncid);
    global_dust_emis_rate_new_param(p) = sum(sum(dust_flux.*area.*year_in_seconds))/10^9; %global dust emission rate in Tg/year
end 

%Calibrating the global dust flux
year_save = year;
%calibrating to global dust flux in ERA-I AERONET-calibrated CLM/CAM runs
old_param_norm_fact = calib_global_dust_flux_old_param/mean(global_dust_emis_rate_old_param(find(year_save==1995):find(year_save==2011)));
new_param_norm_fact = calib_global_dust_flux_new_param/mean(global_dust_emis_rate_new_param(find(year_save==1995):find(year_save==2011)));

%saving the normalization factors, to be used in plot_PNAS_SI_maps
save('CLM_run_norm_factors.mat','old_param_norm_fact','new_param_norm_fact'); %to be read in by plot_SI_maps.m to plot up dust flux changes

global_dust_emis_rate_old_param_calib = old_param_norm_fact*global_dust_emis_rate_old_param;
global_dust_emis_rate_new_param_calib = new_param_norm_fact*global_dust_emis_rate_new_param;

%this reads in the surface temperature averaged over the 25 years
for i = 1:end_year-start_year+1
    year = start_year + i - 1;
    %reading in ECHAM temperature
    ncid = netcdf.open(strcat(model_dir,'ECHAM forcing data/TPQW_',int2str(year-24),'_',int2str(year),'_avg.nc'),'nowrite');
    varid_TBOT = netcdf.inqVarID(ncid,'TBOT');
    TBOT = netcdf.getVar(ncid,varid_TBOT);
    T_surf(i) = sum(sum(temperature_area.*TBOT))/sum(sum(temperature_area));
    netcdf.close(ncid);
end

%calculating the 25 yr averages
for p=1:(end_year-start_year+1)
    global_dust_emis_rate_old_param_25yr_avg(p) = mean(global_dust_emis_rate_old_param_calib(p:p+running_avg_years-1));
    global_dust_emis_rate_new_param_25yr_avg(p) = mean(global_dust_emis_rate_new_param_calib(p:p+running_avg_years-1));
    year_25yr_avg(p) = mean(year_save(p:p+running_avg_years-1));
end
x_plot = 1:size(global_dust_emis_rate_old_param_25yr_avg,2);
figure(2); clf; plot(x_plot, global_dust_emis_rate_old_param_25yr_avg, 'b'); hold on; plot(x_plot, global_dust_emis_rate_new_param_25yr_avg, 'r');
figure(3); clf; plot(x_plot, global_dust_emis_rate_old_param_25yr_avg/global_dust_emis_rate_old_param_25yr_avg(1), 'b'); hold on; plot(x_plot, global_dust_emis_rate_new_param_25yr_avg/global_dust_emis_rate_new_param_25yr_avg(1), 'r');


T_change = T_surf-T_surf(1);
%calculating the relative change in the global dust flux
Q_rel_change_old_param = (global_dust_emis_rate_old_param_25yr_avg - global_dust_emis_rate_old_param_25yr_avg(1))/global_dust_emis_rate_old_param_25yr_avg(1);
Q_rel_change_new_param = (global_dust_emis_rate_new_param_25yr_avg - global_dust_emis_rate_new_param_25yr_avg(1))/global_dust_emis_rate_new_param_25yr_avg(1);

%plotting up results for the old parameterization with original soil moisture
[ver_intercept_Q_old_param, slope_Q_old_param, ver_intercept_Q_old_param_SE, slope_Q_old_param_SE, covariance_Q_old_param] = linearfit (T_change, Q_rel_change_old_param, 0);

%plotting up results for the old parameterization with original soil moisture
[ver_intercept_Q_new_param, slope_Q_new_param, ver_intercept_Q_new_param_SE, slope_Q_new_param_SE, covariance_Q_new_param] = linearfit (T_change, Q_rel_change_new_param, 0);

%calculating kappy and lambda for the old parameterization with new soil moisture
ver_intercept_Q_old_param
kappa_old_param_mean = slope_Q_old_param
kappa_old_param_SE = slope_Q_old_param_SE
kappa_old_param_LC = slope_Q_old_param-2*kappa_old_param_SE
kappa_old_param_UC = slope_Q_old_param+2*kappa_old_param_SE

%calculating kappy and lambda for the new parameterization with new soil moisture
ver_intercept_Q_new_param
kappa_new_param_mean = slope_Q_new_param
kappa_new_param_SE = slope_Q_new_param_SE
kappa_new_param_LC = slope_Q_new_param-2*kappa_new_param_SE
kappa_new_param_UC = slope_Q_new_param+2*kappa_new_param_SE

%saving data, to be read in by calc_kappa_and_feedback_pdf
save('kappa_sim_results.mat','kappa_old_param_mean','kappa_old_param_SE','kappa_old_param_LC','kappa_old_param_UC','kappa_new_param_mean','kappa_new_param_SE','kappa_new_param_LC','kappa_new_param_UC'); 

%writing out the data, to be read in by Origin
fidData=fopen('DustClimFeedback.txt','wt');
fprintf(fidData,'Year  delT  OldParamDustFlux OldParamRelDelQperc  NewParamDustFlux NewParamRelDelQperc  \n'); %FdPar file header
for i=1:end_year-start_year+1   
    fprintf(fidData,'%4.0f  %1.4e %1.4e %1.4e %1.4e %1.4e  \n', year_25yr_avg(i), T_change(i), global_dust_emis_rate_old_param_25yr_avg(i), 100*Q_rel_change_old_param(i), global_dust_emis_rate_new_param_25yr_avg(i), 100*Q_rel_change_new_param(i));
end
fclose(fidData);
