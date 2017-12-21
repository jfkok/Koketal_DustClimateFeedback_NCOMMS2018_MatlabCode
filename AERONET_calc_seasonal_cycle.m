function [calib_fact, corrcoef_lin, RSME, monthly_corrcoef, model_corr_obs_daily] = AERONET_calc_seasonal_cycle (year_start,year_end,min_months,min_days_per_month)

%script that calculates the seasonal cycle in AOD at AERONET stations,
%the result of which is used to calculate the anomaly in AOD for a given
%measurement in AERONET_calc_longterm_trend

fidSeasonalCycle=fopen('AERONET_seasonal_cycle.txt','wt'); %this file contains the calculated seasonal cycle for each station, to be used to eliminate the seasonal cycle from the long-term trend

create_plots = true;

global angstrom_thr

angstrom_thr = 100; %set to large number in order to use all data. this is necessary to investigate a long-term trend (otherwise unclear waht change in AOD means if #of dusty days per month changes)
require_full_year = false; %sets whether to require at least 1 data point for each month to calculate the annual AOD
eliminate_stations = [33, 37, 48]; %eliminates selected stations because not dust-dominated

%read in the filenames of the AERONET data
load_station_filenames_daily;

no_stations = size(station_filenames,1);
for s=1:no_stations %setting the station names
    station_names(s,1:6) = station_filenames(s,33:38);
end %for, setting the station names
%the names of the stations and their coordinates
station_lat = [16.733,    12.200,    36.508,    13.278,    31.922,    13.777,    13.541,    8.320,    30.855,    26.208,    14.394,    14.06,    11.85,    23.717,    43.577,    28.309,    -28.976,    24.907,    22.967,    24.481,    35.517,    22.79,    24.372,    15.345,    31.626,    9.76,    25.495,    25.513,    23.145,    28.473,    13.217,    -25.899,    -43.25,    31.67,    24.87,    28.482,    31.542,    29.503,    36.705, 35.55, 38.283, -16.108, 35.946, 42.623, 13.775, 26.906, 14.709, 38.553, 34.653, 16.864, 22.79, -19.175, 26.513];
station_lon = [337.065,    358.600,    2.881,    354.066,    34.789,    8.99,    2.665,    4.340,    34.782,    50.609,    343.041,   357.55,    356.25,    344.05,    104.419,    343.501,    139.991,    46.397,    54.300,    54.383,    12.632,    5.530,    54.467,    358.521,    351.844,   1.599,    53.146,    56.325,    53.779,    343.753,    12.024,    139.346,    294.691,    352.401,    67.03,    343.679,    74.325,    34.917,    48.507, 8.683, 109.717, 128.749, 104.137, 76.983, 8.984, 75.806, 16.477, 68.858, 358.101, 335.133, 5.530, 15.914, 80.232];
station_source_id = [2, 1, 1, 1, 3, 1, 1, 1, 3, 3, 1, 1, 1, 1, 4, 2, 5, 3, 3, 3, 1, 1, 3, 1, 1, 1, 3, 3, 3, 2, 1, 5, 6, 1, 4, 2, 4, 3, 3, 1, 4, 5, 4, 4, 1, 4, 1, 4, 1, 2, 2, 7, 4];
%1 = North Africa, 2 = North Atlantic, 3 = Middle East, 4 = Rest of Asia, 5= Australia, 6 = South America

%initializing some variables
no_data_matches=zeros(1,no_stations,12); %holds the number of matches between a given model and AERONET station, for each month
station_monthly_AOT=zeros(1,no_stations,12); %holds the number of matches between a given model and AERONET station, for each month
station_monthly_save=zeros(1,no_stations,12*(year_end-year_start+1)); %saving the AOT for each month to calculate the standard deviation
obs_TotAOT = zeros(no_stations,round(365*(year_end-year_start+1)));
station_AOT = zeros(no_stations,6,round(365*(year_end-year_start+1)));

%reading in the AERONET data
for s=1:no_stations %cycling over stations
    temp = read_obs_AOT_daily (station_filenames(s,:));
    station_AOT(s,1:6,1:size(temp,2)) = temp;
end

j=0;
%extracting the model AOD,matching it to the station AOD, and calculating daily, monthly, and total (~annual) comparisons
for r=1:1 %looping over the model runs
    %initializing matrices
    total_AOT_days = zeros(no_stations,1); %keeps track of the number of days that AERONET data is available for each station
    for y=year_start:year_end %looping over the simulated years
        y
        for m=1:12 %looping over the 12 months
            %extracting the AERONET AOD and the model AOD for the same location
            j = j+1;
            for s=1:no_stations %looping over the number of stations
                A=find(station_AOT(s,1,:)==y); B=A(find(station_AOT(s,2,A)==m)); 
                AOT_days=station_AOT(s,3,B); no_AOT_days = size(AOT_days,3);
                obs_TotAOT(s,total_AOT_days(s)+1:total_AOT_days(s)+no_AOT_days)=station_AOT(s,6,B);
                total_AOT_days(s) = total_AOT_days(s) + no_AOT_days;
                if(size(B,1)>min_days_per_month) %in this case there's station data for the given month and year, and this calculates station AOT and model AOT at station location for the same month
                    no_data_matches(r,s,m) = no_data_matches(r,s,m)+1;
                    AOT_dusty_days(1:size(B,1))=station_AOT(s,6,B);
                    station_monthly_AOT(r,s,m) = station_monthly_AOT(r,s,m)+mean(AOT_dusty_days(1:size(B,1)));
                    station_monthly_save(r,s,j) = mean(AOT_dusty_days(1:size(B,1)));
                end %if, in this case there's station data for the given month and year, and this calculates station AOT and model AOT at station location for the same month
            end %for, looping over the number of stations
        end %for, looping over the 12 months
    end %for looping over the simulated years
    
    %determingin which stations meet quality control criteria, calculating the monthly averages for the seasonal cycle, and the monthly and daily correlations
    i=0;
    for s=1:no_stations %looping over the stations
        if (min(eliminate_stations~=s)) %this helps determine whether a station can be used for the seasonal cycle
            if (sum(no_data_matches(r,s,:))>=min_months && sum(no_data_matches(r,s,:)>0)>1 && (min(no_data_matches(r,s,:))>0 || require_full_year == false)) %for a station to be used: (1) must exceed min no of months, (2) must be at least 2 different months, (3) must be at least one data point per month if the program requires a full monthly cycle
                use_station(s) = true;
                i=i+1;
                seasonal_station(i)=s;
            else
                use_station(s) = false;
            end
        else %in this case, the station has been eliminated
            use_station(s) = false;
        end
        station_no_months(s)=sum(no_data_matches(r,s,:));
        avg_monthly_AOT(s,:) = station_monthly_AOT(r,s,:)./no_data_matches(r,s,:);
        clear A; clear B;
        A(1,1:j) = repmat(avg_monthly_AOT(s,:),1,year_end-year_start+1); %extending the (12 point) seasonal cycle to the full data range
        B(1,1:j) = station_monthly_save(r,s,1:j);
        seas_cycle_corr_AOT(s,1:j)=B-A; %the monthly-averaged AOT, corrected for the seasonal cycle        
        if (use_station(s))
            fprintf(fidSeasonalCycle,'%2.0i ',s);
            fprintf(fidSeasonalCycle,station_names(s,:));
            fprintf(fidSeasonalCycle,' %3.2f  %3.2f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  %2.4f  \n',station_lat(s),station_lon(s),avg_monthly_AOT(s,1),avg_monthly_AOT(s,2),avg_monthly_AOT(s,3),avg_monthly_AOT(s,4),avg_monthly_AOT(s,5),avg_monthly_AOT(s,6),avg_monthly_AOT(s,7),avg_monthly_AOT(s,8),avg_monthly_AOT(s,9),avg_monthly_AOT(s,10),avg_monthly_AOT(s,11),avg_monthly_AOT(s,12));
        end
        obs_avg_TotAOT(r,s) = mean(obs_TotAOT(s,1:total_AOT_days(s)));
    end %for, looping over the stations
end %for, looping over the model runs
total_stations = sum(use_station);

if (create_plots)
    create_seasonal_cycle_plots
end

fclose(fidSeasonalCycle);