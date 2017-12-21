fidData_annavg=fopen('AERONET_modeled_annavg_AOD.txt','wt'); %this file contains the modeled and measured AOD at the stations that meet the quality control criteria
fidData_longterm=fopen('AERONET_modeled_longterm_trend_AOD.txt','wt'); %this file contains the modeled and measured longterm trend in AOD at the stations that meet the quality control criteria

%creating the header
fprintf(fidData_annavg,'Year');
for i=1:total_longterm_stations
    string = strcat(', AERONETst',num2str(i),', a1ZenSt',num2str(i),', a1KokSt',num2str(i));
	fprintf(fidData_annavg,string); %FdPar file header
end
fprintf(fidData_annavg,'\n');

%writing out the data
for j=1:year_end-year_start+1 %cycling over all the years
    year=0.5+year_start+j-1; %the year, taking in the middle of the calendar year
    fprintf(fidData_annavg,'%4.1f,  ',year);
    for i=1:total_longterm_stations %cycling over all the stations
        s = has_longterm_trend(i);
        fprintf(fidData_annavg,'%1.4e,  %1.4e,   %1.4e,  ', meas_annual_mean_AOT_seas_corr(s,j),model_annual_mean_AOT_seas_corr(1,s,j),model_annual_mean_AOT_seas_corr(2,s,j));
    end 
    fprintf(fidData_annavg,'\n');
end
fclose(fidData_annavg);

%creating the header
fprintf(fidData_longterm,'Station, NoYears, AODtrendMeasured, AODtrendMeasuredSE, AODtrendA1Zen, AODtrendA1ZenSE, AODtrendA1Kok, AODtrendA1KokSE \n');
for i=1:total_longterm_stations
	s = has_longterm_trend(i);
    string = station_names(s,:);
    fprintf(fidData_longterm,string);
    fprintf(fidData_longterm,', %2.0f,   %1.4e,  %1.4e,  %1.4e,   %1.4e,  %1.4e,   %1.4e \n', noYears(s), 10*AER_slope(s),10*AER_SE_slope(s),10*mod_slope(1,s),10*mod_SE_slope(1,s),10*mod_slope(2,s),10*mod_SE_slope(2,s));
end
