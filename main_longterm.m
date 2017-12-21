function main_longterm

%main script, which calculates the comparison of the CESM AERONET-BASE and
%AERONET-PHYS simulations against AOD measured at AERONET stations

global angstrom_thr

%quality control parameters from code used for 2014 ACP2 paper
min_days_per_month = 5; %minimum dusty days per month for the month to be taken into account; 10 in ACP2 for seasonal cycle comparison
min_months = 4; %the minimum number of total months of data for the station's data to be used; wsa 4 in ACP2
min_months_per_year = 4; %minimum number of months required to take the year into account for the annually-average AOD and the assessment of the long-term AOD trend
angstrom_thr = 100; %set to large number in order to use all data. this is different from ACP2, but necessary to investigate a long-term trend (otherwise unclear waht change in AOD means if #of dusty days per month changes)

%new quality control parameters for assessing long-term trend in data
min_years_separation = 5; %minimum number of years between minimum and maximum year in record for the longterm trend to be used
min_years = 5; %minimum number of years for which there must be data for the longterm trend to be used
min_AOD_angstrom_corr = -sqrt(0.50); %minimum (negative) correlation between AOD and angstrom exponent for the station to be used for long-term trend; set to sqrt(0.50) such that the correlation must explain greater than half of the variance for the station to pass quality control
all_dust_angstrom_thr = 0.4; %threshold of annually-averaged angstrom exponent below which site is considered almost 100% dust and we would not expect a correlation between Angstrom exponent and AOD
dust_AOT_fraction_threshold = 0.75; %minimum fraction of AOD that must be due to dust - as quantified by the average from the four ACP2 simulations - in order for the station to be considered for a longterm AOD trend

recalc_seasonal_cycle = true; %sets whether the seasonal cycle is recalculated, or used from the previous run
reanalysis = 'ERA'; %using ERA-Interim data set

if (reanalysis=='ERA')
    year_start = 1995; %should be 1995, the year in which the runs start
    year_end = 2011; %the year in which the run ends
    model_names = ['NatureComm_rev_run_a1_Zend_src';'NatureComm_rev_run_a1_koknosrc'];
    calib_fact = [1.0711, 1.0572]; %for ERA    
end

if (recalc_seasonal_cycle)
    AERONET_calc_seasonal_cycle (year_start,year_end,min_months,min_days_per_month);
end

[mean_corr, mean_RSME] = AERONET_calc_longterm_trend (model_names, calib_fact, year_start, year_end, min_days_per_month, min_months_per_year, min_years, min_years_separation, min_AOD_angstrom_corr, all_dust_angstrom_thr, dust_AOT_fraction_threshold);
1;