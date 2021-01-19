% generate_elevation.m
%
% This script generates a topography file for the UK region for both the
% GCM, RCM and CPM grids used in UKCP, using the NOAA ETOPO2v2 dataset.
%
% I took the data from the BRIDGE servers (/opt/data/topo_datasets/) and
% official documentation on the data can be found here:
% https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2-2001/readme.txt
%

%% Find which machine is being used
curdir = pwd;
if strcmp(curdir(1:14),'/Users/ak0920/') || strcmp(curdir(1:14),'/Volumes/DataD')
    topodatadir = '/Users/ak0920/Data/IMPRES/OtherData/topo/';
else
    if strcmp(curdir(1:14),'/home/bridge/a')
        topodatadir = '/PreProcessedData/'; 
    end
end

% If the script has been run before, load the pre-processed data:
if exist([topodatadir,'ht_HadUK60.mat'],'file')
    load([topodatadir,'ht_GCM.mat']);
    load([topodatadir,'ht_RCM.mat']);
    load([topodatadir,'ht_CPM.mat']);
    load([topodatadir,'ht_ERA5.mat']);
    load([topodatadir,'ht_HadUK1.mat']);
    load([topodatadir,'ht_HadUK12.mat']);
    load([topodatadir,'ht_HadUK60.mat']);
    
else
    %% Load data
    
    % Elevation data
    ht_raw = ncread([topodatadir,'ETOPO2v2c_f4.nc'],'z',[4800 4100],[1000 500]);
    lon_raw = ncread([topodatadir,'ETOPO2v2c_f4.nc'],'x',[4800],[1000]);
    lat_raw = ncread([topodatadir,'ETOPO2v2c_f4.nc'],'y',[4100],[500]);
    [lat_raw,lon_raw]=meshgrid(lat_raw,lon_raw);
    lat_raw = double(lat_raw); lon_raw = double(lon_raw);
    
    % Lat-long data for all other datasets
    generate_UK_latlon
    
    
    %% Interpolate elevation to lower resolution
    % GCM
    ht_GCM = griddata(lon_raw,lat_raw,ht_raw,long_UK_GCM,lat_UK_GCM);
    
    % RCM
    ht_RCM = griddata(lon_raw,lat_raw,ht_raw,long_UK_RCM,lat_UK_RCM);
    
    % CPM
    ht_CPM = griddata(lon_raw,lat_raw,ht_raw,long_UK_CPM,lat_UK_CPM);
    
    % ERA5
    ht_ERA5 = griddata(lon_raw,lat_raw,ht_raw,double(long_UK_ERA5),double(lat_UK_ERA5));
    
    % HadUK-Grid 1km
    ht_HadUK1 = griddata(lon_raw,lat_raw,ht_raw,long_UK_HadUK1,lat_UK_HadUK1);
    
    % HadUK-Grid 12km
    ht_HadUK12 = griddata(lon_raw,lat_raw,ht_raw,long_UK_HadUK12,lat_UK_HadUK12);
    
    % HadUK-Grid 60km
    ht_HadUK60 = griddata(lon_raw,lat_raw,ht_raw,long_UK_HadUK60,lat_UK_HadUK60);
    
    clear ht_raw lon_raw lat_raw
    
    %% Remove points below sea level
    ht_GCM(ht_GCM<0)=0;
    ht_RCM(ht_RCM<0)=0;
    ht_CPM(ht_CPM<0)=0;
    ht_ERA5(ht_ERA5<0)=0;
    ht_HadUK1(ht_HadUK1<0)=0;
    ht_HadUK12(ht_HadUK12<0)=0;
    ht_HadUK60(ht_HadUK60<0)=0;
    
    
    %% Save the data so not to waste time again in the future
    save('ht_CPM.mat','ht_CPM');
    save('ht_RCM.mat','ht_RCM');
    save('ht_GCM.mat','ht_GCM');
    save('ht_ERA5.mat','ht_ERA5');
    save('ht_HadUK1.mat','ht_HadUK1');
    save('ht_HadUK12.mat','ht_HadUK12');
    save('ht_HadUK60.mat','ht_HadUK60');
    
    
end

clear curdir topodatadir