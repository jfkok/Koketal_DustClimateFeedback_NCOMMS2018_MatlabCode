function [mean_corr, mean_RMSE, corrcoef_lin, monthly_corrcoef, model_corr_obs_daily] = AERONET_calc_longterm_trend (model_names, calib_fact, year_start, year_end, min_days_per_month, min_months_per_year, min_years, min_years_separation, min_AOD_angstrom_corr, all_dust_angstrom_thr, dust_AOT_fraction_threshold)

%script that calculates the long-term AOD trend at dusty AERONET stations,
%and compares those against model simulations

%the location of the model runs
model_dir = '../Data/CESM/Daily-averaged/';
model_text = '_extracted_daily_';
no_model_runs = size(model_names,1);

ref_run_name = model_names(1,:);

%read in the filenames of the AERONET data
load_station_filenames_daily;

no_stations_total = size(station_filenames,1);
for s=1:no_stations_total %setting the station names
    station_names(s,1:6) = station_filenames(s,33:38);
end %for, setting the station names
%the names of the stations and their coordinates
station_lat = [16.733,    12.200,    36.508,    13.278,    31.922,    13.777,    13.541,    8.320,    30.855,    26.208,    14.394,    14.06,    11.85,    23.717,    43.577,    28.309,    -28.976,    24.907,    22.967,    24.481,    35.517,    22.79,    24.372,    15.345,    31.626,    9.76,    25.495,    25.513,    23.145,    28.473,    13.217,    -25.899,    -43.25,    31.67,    24.87,    28.482,    31.542,    29.503,    36.705, 35.55, 38.283, -16.108, 35.946, 42.623, 13.775, 26.906, 14.709, 38.553, 34.653, 16.864, 22.79, -19.175, 26.513];
station_lon = [337.065,    358.600,    2.881,    354.066,    34.789,    8.99,    2.665,    4.340,    34.782,    50.609,    343.041,   357.55,    356.25,    344.05,    104.419,    343.501,    139.991,    46.397,    54.300,    54.383,    12.632,    5.530,    54.467,    358.521,    351.844,   1.599,    53.146,    56.325,    53.779,    343.753,    12.024,    139.346,    294.691,    352.401,    67.03,    343.679,    74.325,    34.917,    48.507, 8.683, 109.717, 128.749, 104.137, 76.983, 8.984, 75.806, 16.477, 68.858, 358.101, 335.133, 5.530, 15.914, 80.232];
station_source_id = [2, 1, 1, 1, 3, 1, 1, 1, 3, 3, 1, 1, 1, 1, 4, 2, 5, 3, 3, 3, 1, 1, 3, 1, 1, 1, 3, 3, 3, 2, 1, 5, 6, 1, 4, 2, 4, 3, 3, 1, 4, 5, 4, 4, 1, 4, 1, 4, 1, 2, 2, 7, 4];
%1 = North Africa, 2 = North Atlantic, 3 = Middle East, 4 = Rest of Asia, 5= Australia, 6 = South America
load('AERONET_station_dust_AOD_fraction.mat'); %loading the fraction of AOD that is dust for each station; from the simulations in ACP2 (written out in main_AERONET_compare_daily in the folder with code for the revised paper)

%reading in each station's seasonal cycle
filename = 'AERONET_seasonal_cycle.txt';
fid = fopen(filename); %open the data file and get the file id
i=0; %keeps track of the number of stations used

while (1)
    A=textscan(fid, ' %u %s %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f', 1); %reading in the data in the format of 7 floating point numbers in a row
    if (size(A{1},1)==0) %in this case the end of the file was reached
        break; 
    end %in this case the end of the file was reached
    i=i+1; %keeps track of the number of stations used
    s = A{1}(1); %the station number is the first column
    array_stations(i)=s;
    AOT_seas_cycle(s,1:12)=cell2mat(A(5:16)); %the seasonal AOT at the station
end
fclose(fid);
no_stations=i;
%'fix no_stations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

%determining the station locations in the model grid box
ncid = netcdf.open('currentall2x2.nc','nowrite');
varid_lon = netcdf.inqVarID(ncid,'lon'); model_lon = netcdf.getVar(ncid,varid_lon);
varid_lat = netcdf.inqVarID(ncid,'lat'); model_lat = netcdf.getVar(ncid,varid_lat);
del_lat = 180/(size(model_lat,1)-1);
del_lon = 360/size(model_lon,1);
for p=1:no_stations
    s=array_stations(p);
    station_lat_no(s)=find(model_lat<station_lat(s)&model_lat+del_lat/2>station_lat(s)|(model_lat>station_lat(s)&model_lat-del_lat/2<station_lat(s)));
    if (station_lon(s)<356.25)
        station_lon_no(s)=find(model_lon<=station_lon(s)&model_lon+del_lon/2>station_lon(s)|(model_lon>station_lon(s)&model_lon-del_lon/2<station_lon(s)));
    elseif (station_lon(s) < 358.75)
        station_lon_no(s)=144;
    else
        station_lon_no(s)=1;
    end
end

%initializing some matrices
obs_TotAOT = zeros(no_stations_total,round(365*(year_end-year_start+1)));
station_AOT = zeros(no_stations_total,6,round(365*(year_end-year_start+1)));
station_angstrom = zeros(no_stations_total,5,round(365*(year_end-year_start+1)));
model_TotAOT_monthly_save = zeros(no_model_runs,no_stations_total,12);
nu_data_matches = zeros(no_model_runs,no_stations_total,12);
no_months = (year_end-year_start+1)*12;
obs_AOT_monthly = zeros(no_stations_total,(year_end-year_start+1)*12);
obs_AOT_monthly_corr = zeros(no_stations_total,(year_end-year_start+1)*12);
obs_angstrom_monthly = zeros(no_stations_total,(year_end-year_start+1)*12);
model_DustAOT_monthly = zeros(no_model_runs, no_stations_total, (year_end-year_start+1)*12);
model_TotAOT_monthly = zeros(no_model_runs, no_stations_total, (year_end-year_start+1)*12);
model_TotAOT_monthly_corr = zeros(no_model_runs, no_stations_total, (year_end-year_start+1)*12);
model_OtherAOT_monthly = zeros(no_model_runs, no_stations_total, (year_end-year_start+1)*12);
model_TotAOT_monthly_save = zeros(no_model_runs, no_stations_total, 12);
model_TotAOT_seas_cycle = zeros(no_model_runs, no_stations_total, 12);

%reading in the AERONET data
for p=1:no_stations %cycling over stations
    s=array_stations(p);
    [temp_AOT,temp_angstrom] = read_obs_AOT_daily (station_filenames(s,:));
    station_AOT(s,1:6,1:size(temp_AOT,2)) = temp_AOT; %info on measured AOT; first row is year, 2nd month, 3rd day, 4th AOT440, 5th AOT675, 6th row AOT550
    station_angstrom(s,1:5,1:size(temp_angstrom,2)) = temp_angstrom; %info on measured angstrom exponent, first row is year, 2nd month, 3rd day, 4th 440-870 nm, 5th 440-675 nm
end

%cd '../../../Total dust flux theory/AERONET/CESM output/Daily/' %to help in opening the filenames, which are otherwise too long                
%extracting the model AOD,matching it to the station AOD, and calculating daily, monthly, and total (~annual) comparisons
for r=1:no_model_runs %looping over the model runs
    %initializing matrices
    DustAOT_daily = zeros(30,size(model_lon,1),size(model_lat,1));
    TotAOT_daily = zeros(30,size(model_lon,1),size(model_lat,1));
    OtherAOT_daily = zeros(30,size(model_lon,1),size(model_lat,1));
    
    total_AOT_days = zeros(no_stations_total,1); %keeps track of the number of days that AERONET data is available for each station
    j=0; %keeps track of number of months since start_year
    model_names(r,:)
    for y=year_start:year_end %looping over the simulated years
        y
        profile on;
        for m=1:12 %looping over the 12 months
            j = j+1; %keeps track of number of months since start_year
            end_day = determine_end_day(m);
            for d=1:end_day
                %extracting the model AOD
                filename = strcat(model_dir,model_names(r,:),'/',model_names(r,:),model_text,num2str(y),'-',num2str(m,'%0.2i'),'-',num2str(d,'%0.2i'),'.nc');
                ncid = netcdf.open(filename,'nowrite');
                varid_TotAOT = netcdf.inqVarID(ncid,'AEROD_v'); TotAOT = netcdf.getVar(ncid,varid_TotAOT);
                varid_DustAOT1 = netcdf.inqVarID(ncid,'ODV_DST01'); DustAOT1 = netcdf.getVar(ncid,varid_DustAOT1);
                varid_DustAOT2 = netcdf.inqVarID(ncid,'ODV_DST02'); DustAOT2 = netcdf.getVar(ncid,varid_DustAOT2);
                varid_DustAOT3 = netcdf.inqVarID(ncid,'ODV_DST03'); DustAOT3 = netcdf.getVar(ncid,varid_DustAOT3);
                varid_DustAOT4 = netcdf.inqVarID(ncid,'ODV_DST04'); DustAOT4 = netcdf.getVar(ncid,varid_DustAOT4);            
                netcdf.close(ncid);
                DustAOT = DustAOT1 + DustAOT2 + DustAOT3 + DustAOT4; %obtaining total dust AOT by summing over the 4 bins
                OtherAOT = TotAOT - DustAOT; %AOT due to other aerosols
                DustAOT=calib_fact(r)*DustAOT; %correcting DustAOT with calibration factor
                TotAOT=OtherAOT+DustAOT; %accounting for correction to DustAOT in TotAOT
                DustAOT_daily(d,:,:) = DustAOT;
                TotAOT_daily(d,:,:) = TotAOT;
                OtherAOT_daily(d,:,:) = OtherAOT;
                close all force; %prevents Matlab from slowing down in opening netCDF files. see http://www.mathworks.com/matlabcentral/answers/13803-matlab-slowing-down-while-reading-netcdf
            end %cycling over the days
            
            %extracting the AERONET AOD and the model AOD for the same location
            for p=1:no_stations %looping over the number of stations
                s=array_stations(p);
                A=find(station_AOT(s,1,:)==y); %finding subset of all measurements taking in year y
                B=A(find(station_AOT(s,2,A)==m)); %finding subset of measurements in array A taking in month m
                AOT_days=squeeze(station_AOT(s,3,B)); %finding the days of month m and year y for which measurements are available
                no_AOT_days = size(AOT_days,1); %the total number of measurements available in month m and year y
                obs_TotAOT(s,total_AOT_days(s)+1:total_AOT_days(s)+no_AOT_days)=station_AOT(s,6,B); %adding the observed AOT for each day in year y and month m to all observed daily AOT for this station s
                obs_Angstrom(s,total_AOT_days(s)+1:total_AOT_days(s)+no_AOT_days)=station_angstrom(s,4,B); %adding the observed AOT for each day in year y and month m to all observed daily AOT for this station s
                total_AOT_days(s) = total_AOT_days(s) + no_AOT_days;
                if(size(B,1)>min_days_per_month) %in this case there's station data for the given month and year, and this calculates station AOT and model AOT at station location for the same month
                    model_DustAOT_monthly(r,s,j)=mean(DustAOT_daily(AOT_days,station_lon_no(s),station_lat_no(s)));
                    model_TotAOT_monthly(r,s,j)=mean(TotAOT_daily(AOT_days,station_lon_no(s),station_lat_no(s)));
                    model_OtherAOT_monthly(r,s,j)=mean(OtherAOT_daily(AOT_days,station_lon_no(s),station_lat_no(s)));
                    model_TotAOT_monthly_save(r,s,m)=model_TotAOT_monthly_save(r,s,m)+model_TotAOT_monthly(r,s,j); %used to calculate the seasonal cycle of AOD, which is then subtracted for the seasonally-adjusted longterm trend
                    nu_data_matches(r,s,m) = nu_data_matches(r,s,m)+1; %keeps track of the number of months for which there is sufficient data
                    %model_TotAOT_monthly_corr(r,s,j)=model_TotAOT_monthly(r,s,j)-AOT_seas_cycle(s,m);
                    if (r==1) %calculating the seasonally-corrected observed AOT at each station
                        obs_AOT_monthly(s,j)=mean(station_AOT(s,6,B));
                        obs_AOT_monthly_corr(s,j)=mean(station_AOT(s,6,B))-AOT_seas_cycle(s,m);
                        obs_angstrom_monthly(s,j)=mean(station_angstrom(s,4,B));
                    end  %if, calculating the seasonally-corrected observed AOT at each station
                else %in this case the AOT for this particular month and station is not defined
                    model_TotAOT_monthly(r,s,j)=NaN;
                    model_TotAOT_monthly_corr(r,s,j)=NaN;
                    obs_AOT_monthly_corr(s,j)=NaN;
                    obs_angstrom_monthly(s,j)=NaN;
                end %if, in this case there's station data for the given month and year, and this calculates station AOT and model AOT at station location for the same month
            end %for, looping over the number of stations
        end %for, looping over the 12 months
        profile viewer
        p = profile('info');
        profile off;
        clear functions; clear mex;
    end %for looping over the simulated years
    for p=1:no_stations
        %correcting the modeling data for its seasonal cycle
        s=array_stations(p);
        model_TotAOT_seas_cycle(r,s,:) = model_TotAOT_monthly_save(r,s,:)./nu_data_matches(r,s,:); %calculating the mean concentration by month
        clear A; clear B;
        temp_array(1:12) = model_TotAOT_seas_cycle(r,s,1:12);        
        A(1,1:no_months) = repmat(temp_array,1,year_end-year_start+1); %extending the (12 point) seasonal cycle to the full data range
        B(1,1:no_months) = model_TotAOT_monthly(r,s,1:no_months);
        model_TotAOT_monthly_corr(r,s,1:no_months)=B-A; %the monthly-averaged concentration, corrected for the seasonal cycle
        %calculating the yearly-averaged measured and modeled AOD based on the seasonally-corrected data
        i=0;
        for y=year_start:year_end %looping over the years to calculate each year's average value
            i=i+1;
            A = 12*(i-1)+1:12*i;
            isfinite_array = A(isfinite(obs_AOT_monthly_corr(s,A)));
            if (size(isfinite_array,2)>min_months_per_year) %applying a quality control for the minimum number of months required per year
                if (r==1)
                    meas_annual_mean_AOT_seas_corr(s,i) = mean(obs_AOT_monthly_corr(s,isfinite_array));
                    meas_annual_mean_angstrom(s,i) = mean(obs_angstrom_monthly(s,isfinite_array));
                end
                model_annual_mean_AOT_seas_corr(r,s,i) = mean(model_TotAOT_monthly_corr(r,s,isfinite_array));
            else
                if (r==1)
                    meas_annual_mean_AOT_seas_corr(s,i) = NaN;
                    meas_annual_mean_angstrom(s,i) = NaN;
                end
                model_annual_mean_AOT_seas_corr(r,s,i) = NaN;
            end
        end        
    end
end %for looping over the model runs

i = 0; j = 0; k = 0; q = 0; %keeps track of number of stations that pass quality control and /or have useful long-term trends
small_angstrom_exp = 0; AOD_correlates_with_angstrom_exp = 0; %initializing the arrays holding the stations that meet the respective quality control criteria

    year_array = 1:(year_end-year_start+1); %used in the fitting
    for p=1:no_stations
        s=array_stations(p);
        use_array = ~isnan(meas_annual_mean_AOT_seas_corr(s,:)); %constructing an array with 0 if year cannot be used, and 1 if year can be used

        %calculating correlation between AOD and angstrom exponent. If this correlation is large and negative, then changes in AOD are driven by changes in dust (or sea salt)
        meas_AOT_corr=meas_annual_mean_AOT_seas_corr(s,use_array); %array with measured AOD; used for determining correlation between measured and modeled AOD and measured AOD and angstrom exponent
        meas_angstrom=meas_annual_mean_angstrom(s,use_array); %array with measured angstrom; used for determining correlation between measured AOD and angstrom exponent
        run_avg_AOT_corr(s) = mean(meas_AOT_corr); %the run-average of the annually-averaged AOT, corrected for the seasonal cycle; this should be within two SE of 0. Can check with 2*std(meas_AOT_corr)/sqrt(size(meas_AOT_corr,2))
        run_avg_angstrom(s) = mean(meas_angstrom); %the run-average of the annually-averaged angstrom exponent
        if (size(meas_AOT_corr,2)>1)
            [corr,pvalue]=corrcoef([meas_AOT_corr', meas_angstrom']); %calculating the correlation between measured AOD and angstrom exponent
            AOD_angstrom_corr(s)=corr(1,2); %the correlation in the annually-averaged AOD and angstrom exponent
            AOD_angstrom_pvalue(s)=pvalue(1,2); 
        else %not enough data to determine the correlation coefficient
            AOD_angstrom_corr(s)=NaN; %the correlation in the annually-averaged AOD and angstrom exponent
            AOD_angstrom_pvalue(s)=NaN; 
        end
        
        %checking whether (i) the station has not been eliminated in the quality control process, (ii) the years with data are spaced at least min_years_separation years apart to make a longterm trend meaningful, (iii) there's enough data (at least min_years) to make a trend statistically meaningful, (iv) most of the AOD is due to dust
        if (max([max(find(use_array==1)),0])-max([0,min(find(use_array==1))])>=min_years_separation && sum(use_array) >= min_years && AERONET_station_dust_AOD_fraction(s) > dust_AOT_fraction_threshold)
            q = q+1;
            passed_quality_control(q)=s; %since it passed the quality control, this station can be used for assessing the longterm trend in AOD
            %checking whether the station is dust-dominated, as indicated by a very low AE. For all these station, dust accounts for > 80% of AOD as quantified by AERONET_station_dust_AOD_fraction
            if (run_avg_angstrom(s) <= all_dust_angstrom_thr) 
                j = j+1;
                small_angstrom_exp(j)=s;
            end
            %checking whether changes in AOD are (negatively) correlated with changes in Angstrom exponent, which indicates that coarse (presumably dust) aerosols drive the AOD change)
            if (AOD_angstrom_corr(s)<=min_AOD_angstrom_corr)
                k = k+1;
                AOD_correlates_with_angstrom_exp(k)=s;
            end                
        end

        if (max(small_angstrom_exp==s)==1||max(AOD_correlates_with_angstrom_exp==s)==1) %AOD long-term trend depends on model's ability to capture dust if AE is very small or if AE correlates strongly with AOD
            i=i+1; %keeps track of number of stations with useful long-term trends
            has_longterm_trend(i)=s; %since it passed the quality control, this station can be used for assessing the longterm trend in AOD
            noYears(s) = sum(~isnan(meas_annual_mean_AOT_seas_corr(s,:))); %number of years the station has data
            %std_Taylor_meas(s) = std(meas_annual_mean_seas_corr(s,use_array));
            
            %calculating the slope of the long-term observed AOD and angstrom trends
            x_fit = year_array(use_array);
            y_fit_AOD = meas_annual_mean_AOT_seas_corr(s,use_array);
            y_fit_angstrom = meas_annual_mean_angstrom(s,use_array);
            [AER_intercept(s), AER_slope(s), AER_SE_intercept(s), AER_SE_slope(s), AER_covariance(s)]=linearfit(x_fit,y_fit_AOD,0);
            [angstrom_intercept(s), angstrom_slope(s), angstrom_SE_intercept(s), angstrom_SE_slope(s), angstrom_covariance(s)]=linearfit(x_fit,y_fit_angstrom,0);
            
            %calculating the correlation of modeled and observed AOD
            for r=1:no_model_runs
                y_fit = squeeze(model_annual_mean_AOT_seas_corr(r,s,use_array))';
                [mod_intercept(r,s), mod_slope(r,s), mod_SE_intercept(r,s), mod_SE_slope(r,s), mod_covariance(r,s)]=linearfit(x_fit,y_fit,0);
                RMSE(r,s) = sqrt(mean((squeeze(model_annual_mean_AOT_seas_corr(r,s,use_array))'-meas_annual_mean_AOT_seas_corr(s,use_array)).^2));
                %std_Taylor_model(r,s) = std(model_annual_mean_seas_corr(r,s,use_array));
                mod_AOT_corr=zeros(1,sum(use_array)); %clearing the array, just to be sure
                mod_AOT_corr(1:sum(use_array))=model_annual_mean_AOT_seas_corr(r,s,use_array); %array with modeled AOD; used for determination correlation between measured and modeled AOD
                [corr,pvalue]=corrcoef([meas_AOT_corr', mod_AOT_corr']); %calculating the correlation between measured and modeled AOD
                trend_corrcoef(r,s)=corr(1,2); %the correlation in the annually-averaged modeled and observed AOD
                trend_pvalue(r,s)=pvalue(1,2); 
            end
        else %in this case, there is not enough data for a useful longterm trend
            trend_corrcoef(r,s)=NaN;
            trend_pvalue(r,s)=NaN;
        end
    end
total_longterm_stations = size(has_longterm_trend,2);

for r=1:no_model_runs %calculating the mean correlation coefficient
    mean_corr(r)=mean(trend_corrcoef(r,has_longterm_trend));
    mean_RMSE(r)=mean(RMSE(r,has_longterm_trend));
    RMSE_longterm(r) = sqrt(sum((mod_slope(r,has_longterm_trend)-AER_slope(has_longterm_trend)).^2))/total_longterm_stations;
    RMSE_longterm_normalized(r) = RMSE_longterm(r)/max(AER_slope(has_longterm_trend))-min(AER_slope(has_longterm_trend));
    [corr,pvalue]=corrcoef([AER_slope(has_longterm_trend)', mod_slope(r,has_longterm_trend)']); 
    corr_longterm_trend(r) = corr(1,2);
end
mean_corr
mean_RMSE
RMSE_longterm_normalized

%plotting the long-term records of AOD at the dusty stations, as well as
%their anomolies after subtracting the seasonal cycle
create_longterm_trend_yearly_plots;

write_AERONET_data