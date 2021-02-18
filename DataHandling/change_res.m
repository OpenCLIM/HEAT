function [datainterp] = change_res(data,resolution)
% [datainterp] = change_res(data,resolution)
%
% Interpolates the reanalysis/obs. data onto the resolution of the UKCP18 
% data for the model evaluation. For obs. data this really refers to the
% 1km HadUK-Grid data (changing this to CPM resolution), as other datasets
% are already on same grid as RCM and GCM data.
%
% Outputs:
%   datainterp = the interpolated data
% 
% Inputs:
%   data = the reanalysis/obs. data
%   resolution = the desired UKCP18 resolution ('RCM', 'GCM' or 'CPM')
%

%% Load the data to begin with
load_xyz

% Figure out what resolution is to be interpolated to
if strcmp(resolution,'RCM')
    lats = lat_UK_RCM;
    lons = long_UK_RCM;
elseif strcmp(resolution,'GCM')
    lats = lat_UK_GCM;
    lons = long_UK_GCM;
elseif strcmp(resolution,'CPM') || strcmp(resolution,'2km')
    lats = lat_UK_CPM;
    lons = long_UK_CPM;
    
end


%% Select the reanalysis/obs. lat-long data
% if length(data(:,1,1,1)) == 121
%     lats_data = lat_UK_ERA5;
%     lons_data = long_UK_ERA5;
% else 
    if length(data(:,1,1,1)) == 900
        lats_data = lat_UK_HadUK1;
        lons_data = long_UK_HadUK1;
%     else 
%         if length(data(:,1,1,1)) == 40
%             lats_data = lat_UK_ERA5(49:88,42:83);
%             lons_data = long_UK_ERA5(49:88,42:83);
%         end
    end
% end


%% Interpolate each time step
% Create empty array to fill
datainterp = zeros(length(lons(:,1)),length(lons(1,:)),length(data(1,1,:)));

% Interpolate each time step
for i = 1:length(data(1,1,:))
    datainterp(:,:,i) = griddata(lons_data,lats_data,data(:,:,i),lons,lats,'nearest');
end

