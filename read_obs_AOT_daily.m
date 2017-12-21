function [AOT,angstrom] = read_obs_AOT_daily (filename)

global angstrom_thr

%read in the data
first_line=textread(filename(:),'%s',1,'delimiter','\t');
header=strread(first_line{1},'%s','delimiter',',');
M=textread(filename(:),'%s','delimiter','\n','headerlines',2);
years=size(M,1)/12;

done=false;
i=0;
ii=0;

%extract the data
while (done==false)
    i=i+1;
    temp_array(1:66)=strread(M{i},'%s','delimiter',',');
    if (max(isnan(str2double(temp_array([7,16]))))==0 && str2double(temp_array(38)) <= angstrom_thr && min(str2double(temp_array([7,16])))>0) %checking that the AOT440 and AOT675 data exist. Otherwise, it's not read in.
        ii=ii+1;
        A(ii,1:66)=temp_array;
    end %if, checking that the AOT440 and AOT675 data exist. Otherwise, it's not read in.
    if (i==size(M,1))
    	done=true;
	end
end
no_data_points = ii;
if (ii==0) %checking whether there is any data
    AOT=zeros(5);
else %in this case there is data
    %extract AOD, first row is year, 2nd month, 3rd day, 4th AOT440, 5th AOT675, 6th row AOT550
    AOT = zeros (6,no_data_points); %initializing AOT
    angstrom = zeros (5,no_data_points); %initializing angstrom
    AOT(4,1:no_data_points) = str2double(A(:,16));
    AOT(5,1:no_data_points) = str2double(A(:,7));
    %extract angstrom exponent, first row is year, 2nd month, 3rd day, 4th 440-870 nm, 5th 440-675 nm
    angstrom(4,1:no_data_points) = str2double(A(:,38)); %angstrom exponent of 440-870 nm, after Dubovik et al. (2002)
    angstrom(5,1:no_data_points) = -log(AOT(5,:)./AOT(4,:))/log(675/440); %angstrom exponent of 440-675 nm, for calculating AOT at 550        
    AOT(6,1:no_data_points) = AOT(5,1:no_data_points).*(550/675).^(-angstrom(5,1:no_data_points));
    date(:,1:10)=char(A(:,1));
    for i=1:no_data_points
        AOT(1,i)=str2double(date(i,7:10));
        AOT(2,i)=str2double(date(i,4:5));
        AOT(3,i)=str2double(date(i,1:2));
        angstrom(1,i)=str2double(date(i,7:10));
        angstrom(2,i)=str2double(date(i,4:5));
        angstrom(3,i)=str2double(date(i,1:2));
    end
end %if