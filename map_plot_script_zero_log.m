function map_plot_script_zero_log(data,axis_lim, sep_ord_magn, no_colors_prefactor, no_labels_per_order_magn, figure_no, colorscheme, message)

% script designed to plot a variable that straddles zero, with increases in
% the variable denoted in red and decreases in blue. The axes will be
% symmetric about zero, and logarithmically spaced.
% The variables are:
% -data: this is the data to be plotted, which should have dimensions equal to size(lon,2)xsize(lat,2). It's assumed that the data will span from 0 to 360 degrees longitude, and this script reformats the data to span from -180 to +180 
% -axis_lim: this is the maximum (absolute) value plotted on the color bar
% -sep_ord_magn: this sets the order of magnitude difference between the plotted value of 1 and what numerical value in the plotted field 1 denotes (for instance, setting sep_ord_magn=2 will set the plotted value of 1 to the data value of 1e-2
% -no_colors_prefactor: this should be an integer that sets the number of colors. Setting it equal to 1 results in a number of colors equal to the number of tick labels minus 1. The two colors in the middle will both be white. 
% -no_labels_per_order_magn: integer that sets the number of labels per order of magnitude plotted
% -colorscheme: string that sets the colorscheme in the cbrewer script
% -message: this is a string that is the title of the plot
% -figure_no: sets the figure number to be used

%coordinates for CESM at 2.5x1.9 degrees
lon = [-180:2.5:180]'; %longitude coordinate. This should span from -180 to +180
lat = [-90:(180/95):90]'; %latitude coordinate. This should span from -90 to +90

%clearing the figure
figure(figure_no); clf; 

%this prepares the data to be plotted on the logarithmic scale on both
%sides of zero
min_value=10^(-sep_ord_magn+1); 
data_log=sign(data).*log10(abs(data)); %creating logarithmic scale for both positive and negative values
data_log=data_log+sep_ord_magn*sign(data); %adding 1 or more orders of magnitude separation between the plotted values and the data's numerical values. This is necessary to create a colorbar with the correct labeling
data_log(abs(data)<min_value)=0; %if data is below the minimum value, which corresponds to 10^(-sep_ord_magn), it will be given the value of zero for plotting purposes
data_log(data<-axis_lim)=-log10(axis_lim)-sep_ord_magn; %data outside of the axis_lim will be shown in the most extreme color on either end
data_log(data>axis_lim)=log10(axis_lim)+sep_ord_magn; %data outside of the axis_lim will be shown in the most extreme color on either end

%this reformats the data to span from -180 to +180
data_plot(1:72,:)=data_log(73:144,:);
data_plot(73:144,:)=data_log(1:72,:);
data_plot(145,:)=data_log(73,:); %adding the -180 degrees longitude data point at the end of the array (column 95) for plotting purposes
data_plot = permute(data_plot,[2,1]);

%this calculates the number of colors to be used
colorbar_labels=(-log10(axis_lim)-sep_ord_magn):1/no_labels_per_order_magn:(log10(axis_lim)+sep_ord_magn); %setting the values on the color bar
for i=1:size(colorbar_labels,2) %printing out the color bar labels
    if (colorbar_labels(i)==0)
        colorbar_string{i} = 0;
    elseif (abs(colorbar_labels(i))-abs(sep_ord_magn)<1)
        colorbar_string{i} = sprintf('%1.1g',sign(colorbar_labels(i))*10^(sign(colorbar_labels(i))*colorbar_labels(i)-sep_ord_magn));
    else
        colorbar_string{i} = sprintf('%2.0f',sign(colorbar_labels(i))*10^(sign(colorbar_labels(i))*colorbar_labels(i)-sep_ord_magn));
    end
end

%this calculates the color map, using the cbrewer script: https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab
noColors=no_colors_prefactor*(size(colorbar_labels,2)-1)-2; %calculating the number of values to be used in the cbrewer script; two white colors will manually be added in the middle 
M(noColors:-1:1,1:3) = cbrewer('div', colorscheme, noColors); %obtaining the color scheme, and reversing is so that negative values are blue and positive values are red
M_temp(1:floor(noColors/2),1:3) = M(1:floor(noColors/2),1:3); %keeping the first half of the colors
M_temp(floor(noColors/2)+1,1:3) = [1,1,1]; %adding the first white in the middle of the color scheme
M_temp(floor(noColors/2)+2,1:3) = [1,1,1]; %adding another white in the middle of the color scheme
M_temp((floor(noColors/2)+3):(noColors+2),1:3) = M((floor(noColors/2)+1):noColors,1:3); %adding the second half of the colors
M=M_temp;
colormap(M);

%this plots the data with the set axis limit
datarange = linspace(-log10(axis_lim)-sep_ord_magn,log10(axis_lim)+sep_ord_magn,100);
contourf(lon,lat,data_plot,datarange,'edgecolor','none');
caxis([-log10(axis_lim)-sep_ord_magn log10(axis_lim)+sep_ord_magn]); %setting the axis limit

hold on
load('topo.mat','topo');
size(topo);
topo2 = topo;
topo2(:,1:180)=topo(:,181:360);
topo2(:,181:360)=topo(:,1:180);
contour([-180:179],-90:89,topo2,[0 0],'k')
axis([-180 180 -90 90])
set(gca,'FontSize',24,'XTick',[-180:45:180],'Ytick',[-90:30:90]) %writing the tick labels
set(gca,'layer','top');
set(gca,'XLim',[-180 180],'YLim',[-90 90],...
    'Box','on',...
    'TickDir','in',...
    'TickLength',[.01 .01],...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'YGrid','on',...
    'XGrid','on',...
    'XColor',[0 0 0],...
    'YColor',[0 0 0],...
    'LineWidth',1);
axis equal
box on
set(gca,'XLim',[-180 180],'YLim',[-90 90]);
pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) .95*pos(3) pos(4)]);
%xlabel('Longitude')
%ylabel('Latitude')
title(message,'FontSize',25)
C = colorbar;
set(C,'YTick',colorbar_labels,'YTickLabel',colorbar_string,'FontSize',25); %this replaces the colorbar values with those from the plotted field (which is logarithmic) with its corresponding values on the linear scale

hold off
