% load_UK_latlon.m
%
% Load the latitude and longitude data for each dataset (e.g. UKCP18, ERA5,
% HadUK-Grid). This can be created from scratch using generate_UK_latlon.m
% but that function does not work on VM or DAFNI.

%% Load data to start with
load('PreProcessedData/long_UK_CPM.mat')
load('PreProcessedData/lat_UK_CPM.mat')
load('PreProcessedData/long_UK_RCM.mat')
load('PreProcessedData/lat_UK_RCM.mat')
load('PreProcessedData/long_UK_GCM.mat')
load('PreProcessedData/lat_UK_GCM.mat')
