%main script, which calculates kappa (fractional change in dust loading per
%degree globally-averaged surface T change) and the resulting global dust
%climate feedback (using the constraints on the dust direct radiative
%effect in Kok et al., Nature Geoscience, 2017)

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Setting parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accurate = true; %set to true for obtaining final results
DRE_mean_Tot_NormPDF = -0.1667; %the mean DRE estimate, in W/m2; from final results for NatGeo paper
DRE_SE_Tot_NormPDF = 0.1876; %the standard deviation of the DRE estimate, in W/m2
DRE_mean_SW_NormPDF = -0.4651; %the mean DRE estimate, in W/m2
DRE_SE_SW_NormPDF = 0.1916; %the standard deviation of the DRE estimate, in W/m2
DRE_mean_LW_NormPDF = 0.2989; %the mean DRE estimate, in W/m2
DRE_SE_LW_NormPDF = 0.0769; %the standard deviation of the DRE estimate, in W/m2
if (accurate)
    delta_DRE = 0.0004; %the spacing between DRE values; needed for the integration to get DCF
    no_kappa_values = 5000; %number of kappa values returned from kernel density
    delta_DCF = 0.0002; %bin spacing for the discrete dust-climate feedback, in W/m2/K
else
    delta_DRE = 0.002; %the spacing between DRE values; needed for the integration to get DCF
    no_kappa_values = 1000; %number of kappa values returned from kernel density
    delta_DCF = 0.001; %bin spacing for the discrete dust-climate feedback, in W/m2/K
end
load('kappa_sim_results.mat'); %reads in kappa_new_param_mean, kappa_new_param_SE, kappa_old_param_mean, and kappa_old_param_SE from calc_global_dust_emission_rate_CLM_runs.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   obtaining the DRE and Kappa pdfs  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%obtaining the DRE pdf from a normal distribution fit to NatGeo results
DRE_Tot_x = [(DRE_mean_Tot_NormPDF-4*DRE_SE_Tot_NormPDF):delta_DRE:(DRE_mean_Tot_NormPDF+4*DRE_SE_Tot_NormPDF)];
DRE_Tot_pdf = normpdf(DRE_Tot_x,DRE_mean_Tot_NormPDF,DRE_SE_Tot_NormPDF);
DRE_SW_x = [(DRE_mean_SW_NormPDF-4*DRE_SE_SW_NormPDF):delta_DRE:(DRE_mean_SW_NormPDF+4*DRE_SE_SW_NormPDF)];
DRE_SW_pdf = normpdf(DRE_SW_x,DRE_mean_SW_NormPDF,DRE_SE_SW_NormPDF);
DRE_LW_x = [(DRE_mean_LW_NormPDF-4*DRE_SE_LW_NormPDF):delta_DRE:(DRE_mean_LW_NormPDF+4*DRE_SE_LW_NormPDF)];
DRE_LW_pdf = normpdf(DRE_LW_x,DRE_mean_LW_NormPDF,DRE_SE_LW_NormPDF);

%obtaining the new parameterization's kappa pdf
delta_kappa_CLM_new = kappa_new_param_SE/12.5;
kappa_new_param_x = [(kappa_new_param_mean-4*kappa_new_param_SE):delta_kappa_CLM_new:(kappa_new_param_mean+4*kappa_new_param_SE)];
kappa_new_param_pdf = normpdf(kappa_new_param_x,kappa_new_param_mean,kappa_new_param_SE);

%obtaining the old parameterization's kappa pdf
delta_kappa_CLM_old = kappa_old_param_SE/12.5;
kappa_old_param_x = [(kappa_old_param_mean-4*kappa_old_param_SE):delta_kappa_CLM_old:(kappa_old_param_mean+4*kappa_old_param_SE)];
kappa_old_param_pdf = normpdf(kappa_old_param_x,kappa_old_param_mean,kappa_old_param_SE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    Calculating kappa and DCF for CMIP5 ensemble  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_models = 18; %number of CMIP5 models with reported change in dust loading
%raw numbers from Supplementary Table 1 in Allen et al. (2016)
CMIP5_RCP85_dust_loading_change = 0.01*[6.9, 23.2, 14.4, 0.76, 0.75, -6.1, 16.2, -2.5, 7.3, 7.1, 7.0, 7.1, 3.5, -2.9, 21.7, -17.4, 2.8, 3.3];
CMIP5_RCP85_mean_surface_T_change = [3.9761198,4.0822472,4.1496222,2.6150008,2.4450878,2.9304459,2.7230816,4.3968033,4.1746206,4.2745423,4.0806194,3.2393152,4.7736550,4.6182339,3.0645104,3.3168087,3.0254152,3.1216775]; %surface temperature changes, provided by Robert Allen by email on 5/1/17
CMIP5_RCP85_kappa = CMIP5_RCP85_dust_loading_change./CMIP5_RCP85_mean_surface_T_change;
%corresponding models: [ACCESS1.0, CanESM2, GFDL-CM3, GFDL-ESM2G,GFDL-ESM2G, GISS-E2-H-p3, GISS-E2-R-p3, HadGEM2-CC, HadGEM2-ES,IPSL-CM5A-LR, IPSL-CM5A-MR, IPSL-CM5B-LR, MIROC-ESM-CHEM, MIROC-ESM,MIROC5, MRI-CGCM3, NorESM1-M, NorESM1-ME, CAM5]

%calculating the pdf of kappa for CMIP5:
A = min([std(CMIP5_RCP85_kappa),iqr(CMIP5_RCP85_kappa)/1.34]); %Eq. (3.30) in Silverman (1986)
n = size(CMIP5_RCP85_kappa,2); %number of models
bw = 0.9*A/(n^(1/5)); %Eq. (3.31) in Silverman (1986)
[kd_kappa,x_kappa,bw_kappa] = ksdensity(CMIP5_RCP85_kappa,'kernel','normal','npoints',no_kappa_values,'bandwidth',bw); x0_kappa = x_kappa(1)-kd_kappa(1)*(x_kappa(2)-x_kappa(1))/(kd_kappa(2)-kd_kappa(1)); 
delta_kappa = (max(x_kappa)-min(x_kappa))/size(x_kappa,2); %the spacing between kappa values; needed for the integration performed to get DCF
figure(2); clf; plot([x0_kappa,x_kappa],[0,kd_kappa],'k-'); hold on;

%calculating the pdf of the dust-climate feedback (DCF) for CMIP5:
k = 0; k_SW = 0; k_LW = 0; %keeps track of all possible values of i and j
for i=1:size(x_kappa,2) %calculates the probability of occurrence for each possible combinatin of discrete DRE and kappa values
    for j=1:size(DRE_Tot_x,2)
        k = k+1;
        DCF_Tot_CMIP5_x(k) = x_kappa(i)*DRE_Tot_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_Tot_CMIP5_prob(k) = (delta_kappa*kd_kappa(i))*(delta_DRE*DRE_Tot_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end
    for j=1:size(DRE_SW_x,2)
        k_SW = k_SW+1;
        DCF_SW_CMIP5_x(k_SW) = x_kappa(i)*DRE_SW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_SW_CMIP5_prob(k_SW) = (delta_kappa*kd_kappa(i))*(delta_DRE*DRE_SW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end
    for j=1:size(DRE_LW_x,2)
        k_LW = k_LW+1;
        DCF_LW_CMIP5_x(k_LW) = x_kappa(i)*DRE_LW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_LW_CMIP5_prob(k_LW) = (delta_kappa*kd_kappa(i))*(delta_DRE*DRE_LW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
end
DCF_Tot_CMIP5_bin_limit = min(DCF_Tot_CMIP5_x):delta_DCF:max(DCF_Tot_CMIP5_x); %bin limits for DCF values, total, for CMIP5
DCF_SW_CMIP5_bin_limit = min(DCF_SW_CMIP5_x):delta_DCF:max(DCF_SW_CMIP5_x); %bin limits for DCF values, SW, for CMIP5
DCF_LW_CMIP5_bin_limit = min(DCF_LW_CMIP5_x):delta_DCF:max(DCF_LW_CMIP5_x); %bin limits for DCF values, LW, for CMIP5
for i=1:size(DCF_Tot_CMIP5_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_Tot_CMIP5_x>DCF_Tot_CMIP5_bin_limit(i)&DCF_Tot_CMIP5_x<DCF_Tot_CMIP5_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_Tot_CMIP5_pdf(i) = sum(DCF_Tot_CMIP5_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
for i=1:size(DCF_SW_CMIP5_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_SW_CMIP5_x>DCF_SW_CMIP5_bin_limit(i)&DCF_SW_CMIP5_x<DCF_SW_CMIP5_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_SW_CMIP5_pdf(i) = sum(DCF_SW_CMIP5_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
for i=1:size(DCF_LW_CMIP5_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_LW_CMIP5_x>DCF_LW_CMIP5_bin_limit(i)&DCF_LW_CMIP5_x<DCF_LW_CMIP5_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_LW_CMIP5_pdf(i) = sum(DCF_LW_CMIP5_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_Tot_CMIP5_binned_x = DCF_Tot_CMIP5_bin_limit(1:end-1)+delta_DCF/2;
delta_DCF_SW_CMIP5_binned_x = DCF_SW_CMIP5_bin_limit(1:end-1)+delta_DCF/2;
delta_DCF_LW_CMIP5_binned_x = DCF_LW_CMIP5_bin_limit(1:end-1)+delta_DCF/2;
figure(3); clf; subplot(3,1,1); plot(delta_DCF_Tot_CMIP5_binned_x,DCF_Tot_CMIP5_pdf); hold on; plot(delta_DCF_SW_CMIP5_binned_x,DCF_SW_CMIP5_pdf,'g'); plot(delta_DCF_LW_CMIP5_binned_x,DCF_LW_CMIP5_pdf,'r'); %plotting the DCF

%calculating the mean and 95% confidence interval of kappa and DCF for the CMIP5 ensemble:
cumpdf_DCF_Tot_CMIP5 = cumsum(DCF_Tot_CMIP5_pdf)/sum(DCF_Tot_CMIP5_pdf);
lower_DCF_CL_CMIP5_Tot = delta_DCF_Tot_CMIP5_binned_x(max(find(cumpdf_DCF_Tot_CMIP5<0.025)));
upper_DCF_CL_CMIP5_Tot = delta_DCF_Tot_CMIP5_binned_x(min(find(cumpdf_DCF_Tot_CMIP5>0.975)));
median_DCF_CMIP5_Tot = delta_DCF_Tot_CMIP5_binned_x(min(find(cumpdf_DCF_Tot_CMIP5>0.50)));
mean_DCF_CMIP5_Tot = sum(delta_DCF_Tot_CMIP5_binned_x.*DCF_Tot_CMIP5_pdf*delta_DCF);
cumpdf_DCF_SW_CMIP5 = cumsum(DCF_SW_CMIP5_pdf)/sum(DCF_SW_CMIP5_pdf);
lower_CL_CMIP5_SW = delta_DCF_SW_CMIP5_binned_x(max(find(cumpdf_DCF_SW_CMIP5<0.025)));
upper_CL_CMIP5_SW = delta_DCF_SW_CMIP5_binned_x(min(find(cumpdf_DCF_SW_CMIP5>0.975)));
mean_DCF_CMIP5_SW = sum(delta_DCF_SW_CMIP5_binned_x.*DCF_SW_CMIP5_pdf*delta_DCF);
cumpdf_DCF_LW_CMIP5 = cumsum(DCF_LW_CMIP5_pdf)/sum(DCF_LW_CMIP5_pdf);
lower_CL_CMIP5_LW = delta_DCF_LW_CMIP5_binned_x(max(find(cumpdf_DCF_LW_CMIP5<0.025)));
upper_CL_CMIP5_LW = delta_DCF_LW_CMIP5_binned_x(min(find(cumpdf_DCF_LW_CMIP5>0.975)));
mean_DCF_CMIP5_LW = sum(delta_DCF_LW_CMIP5_binned_x.*DCF_LW_CMIP5_pdf*delta_DCF);
cumpdf_kappa_CMIP5 = cumsum(kd_kappa)/sum(kd_kappa);
lower_CL_kappa_CMIP5 = x_kappa(max(find(cumpdf_kappa_CMIP5<0.025)));
upper_CL_kappa_CMIP5 = x_kappa(min(find(cumpdf_kappa_CMIP5>0.975)));
median_kappa_CMIP5 = x_kappa(min(find(cumpdf_kappa_CMIP5>0.50)));
mean_kappa_CMIP5 = sum(x_kappa.*kd_kappa*delta_kappa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    Calculating kappa and DCF for CLM model results with old parameterization %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculating the pdf of the dust-climate feedback (DCF) for CLM simulations with the old parameterization:
k = 0; k_SW = 0; k_LW = 0; %keeps track of all possible values of i and j
for i=1:size(kappa_old_param_x,2) %calculates the probability of occurrence for each possible combinatin of discrete DRE and kappa values
    for j=1:size(DRE_Tot_x,2)
        k = k+1;
        DCF_Tot_old_param_x(k) = kappa_old_param_x(i)*DRE_Tot_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_Tot_old_param_prob(k) = (delta_kappa_CLM_old*kappa_old_param_pdf(i))*(delta_DRE*DRE_Tot_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
    for j=1:size(DRE_SW_x,2)
        k_SW = k_SW+1;
        DCF_SW_old_param_x(k_SW) = kappa_old_param_x(i)*DRE_SW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_SW_old_param_prob(k_SW) = (delta_kappa_CLM_old*kappa_old_param_pdf(i))*(delta_DRE*DRE_SW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
    for j=1:size(DRE_LW_x,2)
        k_LW = k_LW+1;
        DCF_LW_old_param_x(k_LW) = kappa_old_param_x(i)*DRE_LW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_LW_old_param_prob(k_LW) = (delta_kappa_CLM_old*kappa_old_param_pdf(i))*(delta_DRE*DRE_LW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
end
DCF_Tot_old_param_bin_limit = min(DCF_Tot_old_param_x):delta_DCF:max(DCF_Tot_old_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_Tot_old_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_Tot_old_param_x>DCF_Tot_old_param_bin_limit(i)&DCF_Tot_old_param_x<DCF_Tot_old_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_Tot_old_param_pdf(i) = sum(DCF_Tot_old_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_Tot_old_param_binned_x = DCF_Tot_old_param_bin_limit(1:end-1)+delta_DCF/2;
DCF_SW_old_param_bin_limit = min(DCF_SW_old_param_x):delta_DCF:max(DCF_SW_old_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_SW_old_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_SW_old_param_x>DCF_SW_old_param_bin_limit(i)&DCF_SW_old_param_x<DCF_SW_old_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_SW_old_param_pdf(i) = sum(DCF_SW_old_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_SW_old_param_binned_x = DCF_SW_old_param_bin_limit(1:end-1)+delta_DCF/2;
DCF_LW_old_param_bin_limit = min(DCF_LW_old_param_x):delta_DCF:max(DCF_LW_old_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_LW_old_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_LW_old_param_x>DCF_LW_old_param_bin_limit(i)&DCF_LW_old_param_x<DCF_LW_old_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_LW_old_param_pdf(i) = sum(DCF_LW_old_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_LW_old_param_binned_x = DCF_LW_old_param_bin_limit(1:end-1)+delta_DCF/2;
figure(3); subplot(3,1,2); plot(delta_DCF_Tot_old_param_binned_x,DCF_Tot_old_param_pdf,'k'); hold on; plot(delta_DCF_SW_old_param_binned_x,DCF_SW_old_param_pdf,'g'); plot(delta_DCF_LW_old_param_binned_x,DCF_LW_old_param_pdf,'r'); %plotting the DCF

%calculating the mean and 95% confidence interval of DCF for CLM results with the old parameterization:
cumpdf_DCF_Tot_old_param = cumsum(DCF_Tot_old_param_pdf)/sum(DCF_Tot_old_param_pdf);
lower_CL_DCF_Tot_old_param = delta_DCF_Tot_old_param_binned_x(max(find(cumpdf_DCF_Tot_old_param<0.025)));
upper_CL_DCF_Tot_old_param = delta_DCF_Tot_old_param_binned_x(min(find(cumpdf_DCF_Tot_old_param>0.975)));
median_DCF_Tot_old_param = delta_DCF_Tot_old_param_binned_x(min(find(cumpdf_DCF_Tot_old_param>0.50)));
mean_DCF_Tot_old_param = sum(delta_DCF_Tot_old_param_binned_x.*DCF_Tot_old_param_pdf*delta_DCF);
cumpdf_DCF_SW_old_param = cumsum(DCF_SW_old_param_pdf)/sum(DCF_SW_old_param_pdf);
lower_CL_DCF_SW_old_param = delta_DCF_SW_old_param_binned_x(max(find(cumpdf_DCF_SW_old_param<0.025)));
upper_CL_DCF_SW_old_param = delta_DCF_SW_old_param_binned_x(min(find(cumpdf_DCF_SW_old_param>0.975)));
median_DCF_SW_old_param = delta_DCF_SW_old_param_binned_x(min(find(cumpdf_DCF_SW_old_param>0.5)));
mean_DCF_SW_old_param = sum(delta_DCF_SW_old_param_binned_x.*DCF_SW_old_param_pdf*delta_DCF);
cumpdf_DCF_LW_old_param = cumsum(DCF_LW_old_param_pdf)/sum(DCF_LW_old_param_pdf);
lower_CL_DCF_LW_old_param = delta_DCF_LW_old_param_binned_x(max(find(cumpdf_DCF_LW_old_param<0.025)));
upper_CL_DCF_LW_old_param = delta_DCF_LW_old_param_binned_x(min(find(cumpdf_DCF_LW_old_param>0.975)));
median_DCF_LW_old_param = delta_DCF_LW_old_param_binned_x(min(find(cumpdf_DCF_LW_old_param>0.5)));
mean_DCF_DCF_LW_old_param = sum(delta_DCF_LW_old_param_binned_x.*DCF_LW_old_param_pdf*delta_DCF);
cumpdf_kappa_old_param = cumsum(kappa_old_param_pdf)/sum(kappa_old_param_pdf);
lower_CL_kappa_old_param = kappa_old_param_x(max(find(cumpdf_kappa_old_param<0.025)));
upper_CL_kappa_old_param = kappa_old_param_x(min(find(cumpdf_kappa_old_param>0.975)));
median_kappa_old_param = kappa_old_param_x(min(find(cumpdf_kappa_old_param>0.50)));
mean_kappa_old_param = sum(kappa_old_param_x.*kappa_old_param_pdf*delta_kappa_CLM_old);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    Calculating kappa and DCF for CLM model results with new parameterization %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculating the pdf of the dust-climate feedback (DCF) for CLM simulations with the new parameterization:
k = 0; k_SW = 0; k_LW = 0; %keeps track of all possible values of i and j
for i=1:size(kappa_new_param_x,2) %calculates the probability of occurrence for each possible combinatin of discrete DRE and kappa values
    for j=1:size(DRE_Tot_x,2)
        k = k+1;
        DCF_Tot_new_param_x(k) = kappa_new_param_x(i)*DRE_Tot_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_Tot_new_param_prob(k) = (delta_kappa_CLM_new*kappa_new_param_pdf(i))*(delta_DRE*DRE_Tot_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
    for j=1:size(DRE_SW_x,2)
        k_SW = k_SW+1;
        DCF_SW_new_param_x(k_SW) = kappa_new_param_x(i)*DRE_SW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_SW_new_param_prob(k_SW) = (delta_kappa_CLM_new*kappa_new_param_pdf(i))*(delta_DRE*DRE_SW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
    for j=1:size(DRE_LW_x,2)
        k_LW = k_LW+1;
        DCF_LW_new_param_x(k_LW) = kappa_new_param_x(i)*DRE_LW_x(j); %the value of the DCF at this particular combination of DRE and kappa values
        DCF_LW_new_param_prob(k_LW) = (delta_kappa_CLM_new*kappa_new_param_pdf(i))*(delta_DRE*DRE_LW_pdf(j)); %the probability of occurrence for this particular combination of DRE and kappa values; sum(DCF_pdf) = 1
    end 
end
DCF_Tot_new_param_bin_limit = min(DCF_Tot_new_param_x):delta_DCF:max(DCF_Tot_new_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_Tot_new_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_Tot_new_param_x>DCF_Tot_new_param_bin_limit(i)&DCF_Tot_new_param_x<DCF_Tot_new_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_Tot_new_param_pdf(i) = sum(DCF_Tot_new_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_Tot_new_param_binned_x = DCF_Tot_new_param_bin_limit(1:end-1)+delta_DCF/2;
DCF_SW_new_param_bin_limit = min(DCF_SW_new_param_x):delta_DCF:max(DCF_SW_new_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_SW_new_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_SW_new_param_x>DCF_SW_new_param_bin_limit(i)&DCF_SW_new_param_x<DCF_SW_new_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_SW_new_param_pdf(i) = sum(DCF_SW_new_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_SW_new_param_binned_x = DCF_SW_new_param_bin_limit(1:end-1)+delta_DCF/2;
DCF_LW_new_param_bin_limit = min(DCF_LW_new_param_x):delta_DCF:max(DCF_LW_new_param_x); %bin limits for DCF values, for CMIP5
for i=1:size(DCF_LW_new_param_bin_limit,2)-1 %this cycles over the DCF bins, and calculates the probability of the DCF being in that particular bin. it does by summing the probabilities of each combination of kappa and DRE values that yields a DCF within the bin
    bin_array = find(DCF_LW_new_param_x>DCF_LW_new_param_bin_limit(i)&DCF_LW_new_param_x<DCF_LW_new_param_bin_limit(i+1)); %finding the array of values with DCF within the bin limits
    DCF_LW_new_param_pdf(i) = sum(DCF_LW_new_param_prob(bin_array))/delta_DCF; %the integral over DCF_CMIP5_pdf must equal 1, so sum(DCF_CMIP5_pdf)*delta_CDF = 1. Since sum(DCF_pdf) = 1, we must have that DCF_CMIP5_pdf(i) = sum(DCF_pdf(bin_array))/delta_DCF
end
delta_DCF_LW_new_param_binned_x = DCF_LW_new_param_bin_limit(1:end-1)+delta_DCF/2;
figure(3); subplot(3,1,3); plot(delta_DCF_Tot_new_param_binned_x,DCF_Tot_new_param_pdf,'k'); hold on; plot(delta_DCF_SW_new_param_binned_x,DCF_SW_new_param_pdf,'g'); plot(delta_DCF_LW_new_param_binned_x,DCF_LW_new_param_pdf,'r'); %plotting the DCF

%calculating the mean and 95% confidence interval of DCF for CLM results with the new parameterization:
cumpdf_DCF_Tot_new_param = cumsum(DCF_Tot_new_param_pdf)/sum(DCF_Tot_new_param_pdf);
lower_CL_DCF_Tot_new_param = delta_DCF_Tot_new_param_binned_x(max(find(cumpdf_DCF_Tot_new_param<0.025)));
upper_CL_DCF_Tot_new_param = delta_DCF_Tot_new_param_binned_x(min(find(cumpdf_DCF_Tot_new_param>0.975)));
median_DCF_Tot_new_param = delta_DCF_Tot_new_param_binned_x(min(find(cumpdf_DCF_Tot_new_param>0.5)));
mean_DCF_Tot_new_param = sum(delta_DCF_Tot_new_param_binned_x.*DCF_Tot_new_param_pdf*delta_DCF);
cumpdf_DCF_SW_new_param = cumsum(DCF_SW_new_param_pdf)/sum(DCF_SW_new_param_pdf);
lower_CL_DCF_SW_new_param = delta_DCF_SW_new_param_binned_x(max(find(cumpdf_DCF_SW_new_param<0.025)));
upper_CL_DCF_SW_new_param = delta_DCF_SW_new_param_binned_x(min(find(cumpdf_DCF_SW_new_param>0.975)));
mean_DCF_SW_new_param = sum(delta_DCF_SW_new_param_binned_x.*DCF_SW_new_param_pdf*delta_DCF);
cumpdf_DCF_LW_new_param = cumsum(DCF_LW_new_param_pdf)/sum(DCF_LW_new_param_pdf);
lower_CL_DCF_LW_new_param = delta_DCF_LW_new_param_binned_x(max(find(cumpdf_DCF_LW_new_param<0.025)));
upper_CL_DCF_LW_new_param = delta_DCF_LW_new_param_binned_x(min(find(cumpdf_DCF_LW_new_param>0.975)));
mean_DCF_LW_new_param = sum(delta_DCF_LW_new_param_binned_x.*DCF_LW_new_param_pdf*delta_DCF);
cumpdf_kappa_new_param = cumsum(kappa_new_param_pdf)/sum(kappa_new_param_pdf);
lower_CL_kappa_new_param = kappa_new_param_x(max(find(cumpdf_kappa_new_param<0.025)));
upper_CL_kappa_new_param = kappa_new_param_x(min(find(cumpdf_kappa_new_param>0.975)));
median_kappa_new_param = kappa_new_param_x(min(find(cumpdf_kappa_new_param>0.50)));
mean_kappa_new_param = sum(kappa_new_param_x.*kappa_new_param_pdf*delta_kappa_CLM_new);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Writing out the kernel densities   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('kappa_pdf.mat','kappa_old_param_x','kappa_old_param_pdf','','kappa_new_param_pdf','x_kappa','kd_kappa','lower_CL_kappa_CMIP5','upper_CL_kappa_CMIP5','median_kappa_CMIP5','mean_kappa_CMIP5','lower_CL_kappa_old_param','upper_CL_kappa_old_param','median_kappa_old_param','mean_kappa_old_param','lower_CL_kappa_new_param','upper_CL_kappa_new_param','median_kappa_new_param','mean_kappa_new_param'); %saving the kappa pdf for use in reginal_dust_DRE
fid = fopen('kappa_pdf_CMIP5.txt','wt');
fprintf(fid,'kappa  prob_dens \n');
fprintf(fid,'%1.4e  %1.4e \n', x0_kappa, 0);
for p = 1:size(kd_kappa,2)
	fprintf(fid,'%1.4e  %1.4e \n', x_kappa(p), kd_kappa(p));
end
fclose(fid);

fid = fopen('kappa_pdf_new_param.txt','wt');
fprintf(fid,'kappa  prob_dens \n');
for p = 1:size(kappa_new_param_x,2)
	fprintf(fid,'%1.4e  %1.4e \n', kappa_new_param_x(p), kappa_new_param_pdf(p));
end
fclose(fid);

fid = fopen('kappa_pdf_old_param.txt','wt');
fprintf(fid,'kappa  prob_dens \n');
for p = 1:size(kappa_old_param_x,2)
	fprintf(fid,'%1.4e  %1.4e \n', kappa_old_param_x(p), kappa_old_param_pdf(p));
end
fclose(fid);

fid = fopen('DCF_pdf_CMIP5.txt','wt');
fprintf(fid,'DCF_tot  prob_dens_tot DCF_SW  prob_dens_SW DCF_LW  prob_dens_LW \n');
for p = 1:max([size(DCF_Tot_CMIP5_pdf,2),size(DCF_SW_CMIP5_pdf,2),size(DCF_LW_CMIP5_pdf,2)])
    if (p <= size(delta_DCF_Tot_CMIP5_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_Tot_CMIP5_binned_x(p), DCF_Tot_CMIP5_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_SW_CMIP5_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_SW_CMIP5_binned_x(p), DCF_SW_CMIP5_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_LW_CMIP5_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e \n', delta_DCF_LW_CMIP5_binned_x(p), DCF_LW_CMIP5_pdf(p));
    else 
        fprintf(fid,'10 -10 \n');
    end        
end
fclose(fid);

fid = fopen('DCF_pdf_new_param.txt','wt');
fprintf(fid,'DCF_tot  prob_dens_tot DCF_SW  prob_dens_SW DCF_LW  prob_dens_LW  \n');
for p = 1:max([size(delta_DCF_Tot_new_param_binned_x,2),size(delta_DCF_SW_new_param_binned_x,2),size(delta_DCF_LW_new_param_binned_x,2)])
    if (p <= size(delta_DCF_Tot_new_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_Tot_new_param_binned_x(p), DCF_Tot_new_param_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_SW_new_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_SW_new_param_binned_x(p), DCF_SW_new_param_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_LW_new_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e \n', delta_DCF_LW_new_param_binned_x(p), DCF_LW_new_param_pdf(p));
    else 
        fprintf(fid,'10 -10 \n');
    end        
end
fclose(fid);

fid = fopen('DCF_pdf_old_param.txt','wt');
fprintf(fid,'DCF_tot  prob_dens_tot DCF_SW  prob_dens_SW DCF_LW  prob_dens_LW  \n');
for p = 1:max([size(delta_DCF_Tot_old_param_binned_x,2),size(delta_DCF_SW_old_param_binned_x,2),size(delta_DCF_LW_old_param_binned_x,2)])
    if (p <= size(delta_DCF_Tot_old_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_Tot_old_param_binned_x(p), DCF_Tot_old_param_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_SW_old_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e ', delta_DCF_SW_old_param_binned_x(p), DCF_SW_old_param_pdf(p));
    else 
        fprintf(fid,'10 -10 ');
    end
    if (p <= size(delta_DCF_LW_old_param_binned_x,2))
        fprintf(fid,'%1.4e  %1.4e \n', delta_DCF_LW_old_param_binned_x(p), DCF_LW_old_param_pdf(p));
    else 
        fprintf(fid,'10 -10 \n');
    end        
end
fclose(fid);