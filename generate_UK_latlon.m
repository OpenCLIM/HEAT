% generate_UK_latlon.m
%
% Load the latitude and longitude data for each dataset (e.g. UKCP18, ERA5,
% HadUK-Grid) and merge if necessary so it is a 2D grid rather than a 1D
% vector.
% 
% The RCM and GCM grids are the British National grid so the lats and longs 
% need to be transformed before the areas calculated.

%% Find which machine is being used
curdir = pwd;
if strcmp(curdir(1:14),'/Users/ak0920/') || strcmp(curdir(1:14),'/Volumes/DataD')
    HadUKdir = '/Users/ak0920/Data/IMPRES/OtherData/HadUK-Grid/tas/';
    ncdatadir1 = '/Users/ak0920/Data/IMPRES/';
    ncdatadir2 = '/huss/run01/';
else
    if strcmp(curdir(1:14),'/home/bridge/a')
        HadUKdir = '/export/anthropocene/array-01/ak0920/HadUKGrid/tas';
        ncdatadir1 = '/export/anthropocene/array-01/ak0920/ukcp18_data/';
        ncdatadir2 = '/huss/day/run01/';
    end
end


%% Load UKCP18 CPM, RCM and GCM grid data
% CPM
long_UK_CPM = ncread([ncdatadir1,'2km',ncdatadir2,'huss_rcp85_land-cpm_uk_2.2km_01_day_19801201-19811130.nc'],'longitude');
lat_UK_CPM = ncread([ncdatadir1,'2km',ncdatadir2,'huss_rcp85_land-cpm_uk_2.2km_01_day_19801201-19811130.nc'],'latitude');
% RCM
x_UK_RCM = ncread([ncdatadir1,'12km',ncdatadir2,'huss_rcp85_land-rcm_uk_12km_01_day_19801201-19901130.nc'],'projection_x_coordinate');
y_UK_RCM = ncread([ncdatadir1,'12km',ncdatadir2,'huss_rcp85_land-rcm_uk_12km_01_day_19801201-19901130.nc'],'projection_y_coordinate');
% GCM
x_UK_GCM = ncread([ncdatadir1,'60km',ncdatadir2,'huss_rcp85_land-gcm_uk_60km_01_day_19791201-19891130.nc'],'projection_x_coordinate');
y_UK_GCM = ncread([ncdatadir1,'60km',ncdatadir2,'huss_rcp85_land-gcm_uk_60km_01_day_19791201-19891130.nc'],'projection_y_coordinate');

% Merge the RCM and GCM grids into 2D and convert from OS to lat-long proj.
[xs_RCM,ys_RCM]=meshgrid(x_UK_RCM,y_UK_RCM);
xs_RCM = flipud(rot90(xs_RCM)); 
ys_RCM = rot90(ys_RCM);

% This line should work, but for some reason just repeats the first row
[lat_UK_RCM,long_UK_RCM] = os2ll(xs_RCM,ys_RCM);
% This is a work around:
for i = 1:length(xs_RCM(:,1))
    [lat_UK_RCM(i,:),long_UK_RCM(i,:)] = os2ll(xs_RCM(i,:),ys_RCM(i,:));
end

[xs_GCM,ys_GCM]=meshgrid(x_UK_GCM,y_UK_GCM);
xs_GCM = flipud(rot90(xs_GCM)); 
ys_GCM = rot90(ys_GCM);
[lat_UK_GCM,long_UK_GCM] = os2ll(xs_GCM,ys_GCM);


%% Load ERA5 grid data
lat_UK_ERA5 = double(load_ERA5_data('latitude'));
long_UK_ERA5 = double(load_ERA5_data('longitude'));

% Merge into 2D
[lat_UK_ERA5,long_UK_ERA5]=meshgrid(lat_UK_ERA5,long_UK_ERA5);
lat_UK_ERA5 = fliplr(lat_UK_ERA5);


%% Load HadUK-Grid grid data
% 1km is most important, but 12km (and possibly 60km grid) doesn't have
% quite the same limits at the UKCP18 RCM (and GCM) simulations, so load
% them too:
long_UK_HadUK1 = ncread([HadUKdir,'1km/tas_hadukgrid_uk_1km_mon_198101-198112.nc'],'longitude');
lat_UK_HadUK1 = ncread([HadUKdir,'1km/tas_hadukgrid_uk_1km_mon_198101-198112.nc'],'latitude');

long_UK_HadUK12 = ncread([HadUKdir,'12km/tas_hadukgrid_uk_12km_mon_198101-198112.nc'],'longitude');
lat_UK_HadUK12 = ncread([HadUKdir,'12km/tas_hadukgrid_uk_12km_mon_198101-198112.nc'],'latitude');

long_UK_HadUK60 = ncread([HadUKdir,'60km/tas_hadukgrid_uk_60km_mon_188401-188412.nc'],'longitude');
lat_UK_HadUK60 = ncread([HadUKdir,'60km/tas_hadukgrid_uk_60km_mon_188401-188412.nc'],'latitude');

% Testing difference between method of converting BNG to lat-long
x_UK_HadUK12 = ncread([HadUKdir,'12km/tas_hadukgrid_uk_12km_mon_198101-198112.nc'],'projection_x_coordinate');
y_UK_HadUK12 = ncread([HadUKdir,'12km/tas_hadukgrid_uk_12km_mon_198101-198112.nc'],'projection_y_coordinate');

x_UK_HadUK60 = ncread([HadUKdir,'60km/tas_hadukgrid_uk_60km_mon_188401-188412.nc'],'projection_x_coordinate');
y_UK_HadUK60 = ncread([HadUKdir,'60km/tas_hadukgrid_uk_60km_mon_188401-188412.nc'],'projection_y_coordinate');

% Merge the OSNG grids into 2D and convert from OSNG to lat-long proj.
[xs_HadUK12,ys_HadUK12]=meshgrid(x_UK_HadUK12,y_UK_HadUK12);
xs_HadUK12 = flipud(rot90(xs_HadUK12)); 
ys_HadUK12 = rot90(ys_HadUK12);
[lat_UK_HadUK12a,long_UK_HadUK12a] = os2ll(xs_HadUK12,ys_HadUK12);

[xs_HadUK60,ys_HadUK60]=meshgrid(x_UK_HadUK60,y_UK_HadUK60);
xs_HadUK60 = flipud(rot90(xs_HadUK60)); 
ys_HadUK60 = rot90(ys_HadUK60);
[lat_UK_HadUK60a,long_UK_HadUK60a] = os2ll(xs_HadUK60,ys_HadUK60);


% %% Check everything is orientated correctly
% disp('All data should start in the southwest of the domain:')
% disp(['GCM: ',num2str(lat_UK_GCM(1,1)),'°N, ',num2str(long_UK_GCM(1,1)),'°W'])
% disp(['RCM: ',num2str(lat_UK_RCM(1,1)),'°N, ',num2str(long_UK_RCM(1,1)),'°W'])
% disp(['CPM: ',num2str(lat_UK_CPM(1,1)),'°N, ',num2str(long_UK_CPM(1,1)),'°W'])
% disp(['HadUK1: ',num2str(lat_UK_HadUK1(1,1)),'°N, ',num2str(long_UK_HadUK1(1,1)),'°W'])
% disp(['HadUK12: ',num2str(lat_UK_HadUK12(1,1)),'°N, ',num2str(long_UK_HadUK12(1,1)),'°W'])
% disp(['HadUK60: ',num2str(lat_UK_HadUK60(1,1)),'°N, ',num2str(long_UK_HadUK60(1,1)),'°W'])
% disp(['ERA5: ',num2str(lat_UK_ERA5(1,1)),'°N, ',num2str(long_UK_ERA5(1,1)),'°W'])

%% Tidy up
clear x_UK_RCM y_UK_RCM x_UK_GCM y_UK_GCM x_UK_HadUK12 y_UK_HadUK12 x_UK_HadUK60 y_UK_HadUK60 xs ys
clear curdir HadUKdir ncdatadir1 ncdatadir2
