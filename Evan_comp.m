%this script calculates the comparison against the AVHRR data set of Evan
%and Mukhopadhyay (2010) and Evan et al. (2014)

clear all;

model_dir = '../Data/CESM/AVHRR comparison/'; %location of the data
start_year = 1982;

%reading in the AVHRR data from Amato's netcdf file
file_name = strcat(model_dir,'tnatl.dust.nc'); %file from Amato Evan containing the AVHRR data
ncid = netcdf.open(file_name,'nowrite');
varid_tau_d = netcdf.inqVarID(ncid,'dust_optical_depth'); tau_d_AVHRR_monthly = netcdf.getVar(ncid,varid_tau_d); %yearly averaged dust AOD
varid_time = netcdf.inqVarID(ncid,'time'); time_AVHRR = netcdf.getVar(ncid,varid_time); %yearly averaged dust AOD
varid_lat = netcdf.inqVarID(ncid,'latitude'); lat_AVHRR = netcdf.getVar(ncid,varid_lat); %yearly averaged dust AOD
varid_lon = netcdf.inqVarID(ncid,'longitude'); lon_AVHRR = netcdf.getVar(ncid,varid_lon); %yearly averaged dust AOD
array_lat0 = 1:10;
array_lat1 = 11:20; %corresponds to 10-20 N, so comparable to Evan et al. (2014) results
array_lat2 = 21:30;
array_lon = 36:45;
A0=cos(3.14159*lat_AVHRR(array_lat0)/180);
rel_area0=repmat(A0,1,10);
A1=cos(3.14159*lat_AVHRR(array_lat1)/180);
rel_area1=repmat(A1,1,10);
A2=cos(3.14159*lat_AVHRR(array_lat2)/180);
rel_area2=repmat(A2,1,10);
netcdf.close(ncid);

i=0;
end_year_AVHRR = 2009;
AVHRR_years = (1982+0.5):(end_year_AVHRR+0.5);
AVHRR_tau_d_alt = [0.34389, 0.46997,0.49471,0.50008,0.40082,0.42884,0.41483,0.40678,0.38204,0.42347,0.38264,0.34061,0.35581,0.33674,0.33107,0.35999,0.38651,0.36475,0.36446,0.35462,0.35492,0.31617,0.29709,0.28398,0.33674,0.34478,0.33346];
for y=1982:end_year_AVHRR
    i=i+1;
    tau_d_AVHRR_yearly(:,:) = mean(tau_d_AVHRR_monthly(:,:,(1+12*(y-1982)):12*(y+1-1982)),3);
    AVHRR_tau_d0(i) = mean(mean(rel_area0.*tau_d_AVHRR_yearly(array_lat0,array_lon)));
    AVHRR_tau_d1(i) = mean(mean(rel_area1.*tau_d_AVHRR_yearly(array_lat1,array_lon)));
    AVHRR_tau_d2(i) = mean(mean(rel_area2.*tau_d_AVHRR_yearly(array_lat2,array_lon)));
    Cape_Verde_yearly(i) = mean(tau_d_AVHRR_monthly(17,43,(1+12*(y-1982)):12*(y+1-1982)),3);
end
figure(1); clf; plot(AVHRR_years,AVHRR_tau_d1,'k'); hold on; plot(AVHRR_years,AVHRR_tau_d2,'b'); plot(AVHRR_years,AVHRR_tau_d0,'m'); plot(AVHRR_years(1:end-1),AVHRR_tau_d_alt,'r'); plot(AVHRR_years,Cape_Verde_yearly,'g');

model_names = ['NatureComm_rev_run_a1_kokMerra';'NatureComm_rev_run_a1_ZenMERRA']; 
calib_fact = [0.9561,0.8786]; %calibration factors from tuning against AERONET AOD, per Kok et al. (2014b);
    
nu_model_runs = size(model_names,1);
end_year = 2011;
no_years = end_year - start_year + 1;

AVHRR_start = start_year - AVHRR_years(1) + 1;
AVHRR_end = 27;
CESM_start = 1;
CESM_end = 2008 - start_year + 1;

for r=1:nu_model_runs
    filename = strcat(model_dir,model_names(r,:),'Evan_comp.nc');
    ncid = netcdf.open(filename,'nowrite');
    varid_tau_d_yearly = netcdf.inqVarID(ncid,'total_dust_aod_yearly'); tau_d_yearly = netcdf.getVar(ncid,varid_tau_d_yearly); %yearly averaged dust AOD
    varid_DMP_yearly = netcdf.inqVarID(ncid,'dust_mass_path_yearly'); DMP_yearly = netcdf.getVar(ncid,varid_DMP_yearly); %yearly averaged dust mass path
    varid_Evan_approx_tau_d = netcdf.inqVarID(ncid,'Evan_avg_dust_aod_yearly'); Evan_approx_tau_d = netcdf.getVar(ncid,varid_Evan_approx_tau_d); %yearly averaged dust AOD over approximately the Evan domain
    varid_Evan_approx_DMP = netcdf.inqVarID(ncid,'Evan_avg_dust_mass_path_yearly'); Evan_approx_DMP = netcdf.getVar(ncid,varid_Evan_approx_DMP); %yearly averaged dust AOD over approximately the Evan domain
    netcdf.close(ncid);
    for y=1:no_years
        %calculating the average of tau and DMP over the study domain, accounting for partial grid boxes
        clear temp;
        temp = sum(sum(tau_d_yearly(134:136,55:58,y))); %12 full grid boxes
        temp = temp + 0.5*sum(tau_d_yearly(133,55:58,y)) + 0.5*sum(tau_d_yearly(137,55:58,y)); %8 half grid boxes at W and E boundaries
        temp = temp + 0.5556*sum(tau_d_yearly(134:136,59,y)); %3x 55.56% grid boxes at N boundary
        temp = temp + 0.2778*tau_d_yearly(133,59,y) + 0.2778*tau_d_yearly(137,59,y); %2x 27.78% grid boxes at NW and NE corners
        temp = temp + 0.7223*sum(tau_d_yearly(134:136,54,y)); %3x 77.23% grid boxes at S boundary
        temp = temp + 0.3611*tau_d_yearly(133,54,y) + 0.3611*tau_d_yearly(137,54,y); %2x 36.11% grid boxes at NW and NE corners
        Evan_tau_d1(r,y) = calib_fact(r)*temp/(12+8*0.5+3*0.5556+2*0.2778+3*0.7723+2*0.3611); %the average tau_d in the Evan et al. (GRL, 2014) of 10-20N and 20-30W
        %Evan_tau_d1(r,y) = tau_d_yearly(136,56,y); 'Fix Evan_tau_d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        clear temp;
        temp = sum(sum(DMP_yearly(134:136,55:58,y))); %12 full grid boxes
        temp = temp + 0.5*sum(DMP_yearly(133,55:58,y)) + 0.5*sum(DMP_yearly(137,55:58,y)); %8 half grid boxes at W and E boundaries
        temp = temp + 0.5556*sum(DMP_yearly(134:136,59,y)); %3x 55.56% grid boxes at N boundary
        temp = temp + 0.2778*DMP_yearly(133,59,y) + 0.2778*DMP_yearly(137,59,y); %2x 27.78% grid boxes at NW and NE corners
        temp = temp + 0.7223*sum(DMP_yearly(134:136,54,y)); %3x 77.23% grid boxes at S boundary
        temp = temp + 0.3611*DMP_yearly(133,54,y) + 0.3611*DMP_yearly(137,54,y); %2x 36.11% grid boxes at NW and NE corners
        Evan_DMP1(r,y) = calib_fact(r)*temp/(12+8*0.5+3*0.5556+2*0.2778+3*0.7723+2*0.3611); %the average tau_d in the Evan et al. (GRL, 2014) of 10-20N and 20-30W        
    end
    A_corr = AVHRR_tau_d1(AVHRR_start:AVHRR_end);
    B_corr = Evan_tau_d1(r,CESM_start:CESM_end);
    [corr,pvalue]=corrcoef([A_corr', B_corr']);
    AVHRR_corr(r) = corr(1,2);
    AVHRR_pvalue(r) = pvalue(1,2);
    DMP_start = 2001 - start_year + 1; %starting in 2001
    DMP_end = 2011 - start_year + 1; %ending in 2011
    avg_Evan_DMP1(r) = mean(Evan_DMP1(r,DMP_start:DMP_end));
    avg_Evan_tau_d1(r) = mean(Evan_tau_d1(r,DMP_start:DMP_end));
end
plot_years = start_year:2008;
figure(1); clf; plot(plot_years,AVHRR_tau_d1(AVHRR_start:AVHRR_end),'ko-');

[AVHRR_intercept1, AVHRR_slope1, AVHRR_SE_intercept1, AVHRR_SE_slope1, AVHRR_covariance1]=linearfit(AVHRR_start:AVHRR_end,AVHRR_tau_d1(AVHRR_start:AVHRR_end),0);
AVHRR_mean_tau_d1 = mean(AVHRR_tau_d1(1,AVHRR_start:AVHRR_end)); %mean of the AOD in the study region for the AVHRR measurements

figure(1); hold on; plot(plot_years,AVHRR_intercept1+AVHRR_slope1*(plot_years-start_year),'k--');

figure(1); hold on; plot(plot_years,Evan_tau_d1(1,CESM_start:CESM_end),'go-'); %a=1, Kok
[kok_intercept1, kok_slope1, kok_SE_intercept1, kok_SE_slope1, kok_covariance1]=linearfit(AVHRR_start:AVHRR_end,Evan_tau_d1(1,AVHRR_start:AVHRR_end),0);
Kok_mean_tau_d1 = mean(Evan_tau_d1(1,AVHRR_start:AVHRR_end)); %mean of the AOD in the study region for the Kok simulation
figure(1); hold on; plot(plot_years,kok_intercept1+kok_slope1*(plot_years-start_year),'g--');

figure(1); hold on; plot(plot_years,Evan_tau_d1(2,CESM_start:CESM_end),'bo-'); %a=1, Zender
[Zen_intercept1, Zen_slope1, Zen_SE_intercept1, Zen_SE_slope1, Zen_covariance1]=linearfit(AVHRR_start:AVHRR_end,Evan_tau_d1(2,AVHRR_start:AVHRR_end),0);
Zen_mean_tau_d1 = mean(Evan_tau_d1(2,AVHRR_start:AVHRR_end)); %mean of the AOD in the study region for the Zender simulation
figure(1); hold on; plot(plot_years,Zen_intercept1+Zen_slope1*(plot_years-start_year),'b--'); 
    
fidData=fopen('Evan_output.txt','wt');
fprintf(fidData,'Year  AVHRR   a1_Zender   Kok \n')
for i=1:size(plot_years,2)
	fprintf(fidData,'%4.0i  %1.5f  %1.5f  %1.5f \n', plot_years(i), AVHRR_tau_d1(i), Evan_tau_d1(2,i), Evan_tau_d1(1,i));
end
fclose(fidData);
    
fidFits=fopen('Evan_fits.txt','wt');
fprintf(fidData,'ModelOrObs  Intercept   Slope   Mean    Intercept_SE   Slope_SE    Imag_covariance \n');
fprintf(fidData,'AVHRR1    %1.7f   %1.7f   %1.7f   %1.7f   %1.7f   %1.7f \n', AVHRR_intercept1, AVHRR_slope1, AVHRR_mean_tau_d1, AVHRR_SE_intercept1, AVHRR_SE_slope1, imag(AVHRR_covariance1));
fprintf(fidData,'a1_Zender1    %1.7f   %1.7f   %1.7f   %1.7f   %1.7f   %1.7f \n', Zen_intercept1, Zen_slope1, Zen_mean_tau_d1, Zen_SE_intercept1, Zen_SE_slope1, imag(Zen_covariance1));
fprintf(fidData,'Kok1    %1.7f   %1.7f   %1.7f   %1.7f   %1.7f   %1.7f \n', kok_intercept1, kok_slope1, Kok_mean_tau_d1, kok_SE_intercept1, kok_SE_slope1, imag(kok_covariance1));    
fclose(fidFits);

[AVHRR_slope1, AVHRR_SE_slope1]
[kok_slope1, kok_SE_slope1]
[Zen_slope1, Zen_SE_slope1]