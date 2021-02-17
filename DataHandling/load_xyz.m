% load_xyz.m
%
% Load .mat files providing elevation (for pressure calculations) and
% lat-long area data for spatial averaging. These .mat files are stored in
% PreProcessedData/ and are generate by scripts in Development/ which are
% not included in the main Git package.
%

% Check if the data exists
if ~exist('PreProcessedData','dir')
    disp('Error: PreProcessedData not found ? please contact alan.kennedy@bristol.ac.uk');
    return
end

% Load the elevation data (from generate_elevation.m)
load('PreProcessedData/ht_60km.mat');
load('PreProcessedData/ht_12km.mat');
load('PreProcessedData/ht_2km.mat');
load('PreProcessedData/ht_025deg.mat');
load('PreProcessedData/ht_1km.mat');

% Load the data masks for the UK (from generate_UK_datamask.m)
load('PreProcessedData/datamask2.mat')
load('PreProcessedData/datamask025deg.mat')
load('PreProcessedData/datamask1.mat')
load('PreProcessedData/datamask12.mat')
load('PreProcessedData/datamask60.mat')
datamask025dega = datamask025deg(49:88,42:83);

% Load LSM data
load('PreProcessedData/LSM2.mat')
load('PreProcessedData/LSM1.mat')
load('PreProcessedData/LSM12.mat')
load('PreProcessedData/LSM60.mat')

% Load the lat-long data directly from the data
generate_UK_latlon

% Load the area of each lat-lon box (from generate_UK_latlon_area.m)
load PreProcessedData/areas_60km_abs.mat
load PreProcessedData/areas_60km_frac.mat
load PreProcessedData/areas_60km_frac_UK.mat
load PreProcessedData/areas_12km_abs.mat
load PreProcessedData/areas_12km_frac.mat
load PreProcessedData/areas_12km_frac_UK.mat
load PreProcessedData/areas_2km_abs.mat
load PreProcessedData/areas_2km_frac.mat
load PreProcessedData/areas_2km_frac_UK.mat
load PreProcessedData/areas_1km_abs.mat
load PreProcessedData/areas_1km_frac.mat
load PreProcessedData/areas_1km_frac_UK.mat
load PreProcessedData/areas_025deg_abs.mat
load PreProcessedData/areas_025deg_frac.mat
load PreProcessedData/areas_025deg_frac_UK.mat

