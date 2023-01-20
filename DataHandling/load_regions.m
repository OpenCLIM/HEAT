% load_regions.m
% 
% Load the preprocessed regional masks then create the absolute and
% fractional areas for doing spatial averages etc. 
% 
% The regions have been previously calculated using load_UK_regions.m in
% Development/ (this may require directory structure to be updated).
% 

%% Load data to start with
load('PreProcessedData/UKregions025deg.mat')
load('PreProcessedData/UKregions2.mat')
load('PreProcessedData/UKregions12.mat')
load('PreProcessedData/UKregions1.mat')
load('PreProcessedData/UKregions60.mat')

% Load lat-long area data for whole of UK if required
if ~exist('areas_12km_abs','var')
    load_xyz
end

% Regions for reference:
regs = {'Scotland','North East','North West','Yorkshire and the Humber','East Midlands','West Midlands','East of England','Greater London','South East','South West','Wales','Northern Ireland','Isle of Man'};

%% Calculate the absolute and fractional areas for each region
% Create blank array to fill
areas_025deg_abs_regions = zeros(length(areas_025deg_abs(:,1)),length(areas_025deg_abs(1,:)),12);
areas_025deg_frac_regions = zeros(length(areas_025deg_abs(:,1)),length(areas_025deg_abs(1,:)),12);
areas_2km_abs_regions = zeros(length(areas_2km_abs(:,1)),length(areas_2km_abs(1,:)),12);
areas_2km_frac_regions = zeros(length(areas_2km_abs(:,1)),length(areas_2km_abs(1,:)),12);
areas_12km_abs_regions = zeros(length(areas_12km_abs(:,1)),length(areas_12km_abs(1,:)),12);
areas_12km_frac_regions = zeros(length(areas_12km_abs(:,1)),length(areas_12km_abs(1,:)),12);
areas_60km_abs_regions = zeros(length(areas_60km_abs(:,1)),length(areas_60km_abs(1,:)),12);
areas_60km_frac_regions = zeros(length(areas_60km_abs(:,1)),length(areas_60km_abs(1,:)),12);

% Go through each region
for r=1:12
    areas_025deg_abs_regions(:,:,r) = areas_025deg_abs.*(UKregions025deg == r);
    areas_025deg_frac_regions(:,:,r) = areas_025deg_abs_regions(:,:,r)./nansum(nansum(areas_025deg_abs_regions(:,:,r)));
    areas_2km_abs_regions(:,:,r) = areas_2km_abs.*(UKregions2 == r);
    areas_2km_frac_regions(:,:,r) = areas_2km_abs_regions(:,:,r)./nansum(nansum(areas_2km_abs_regions(:,:,r)));
    areas_12km_abs_regions(:,:,r) = areas_12km_abs.*(UKregions12 == r);
    areas_12km_frac_regions(:,:,r) = areas_12km_abs_regions(:,:,r)./nansum(nansum(areas_12km_abs_regions(:,:,r)));
    areas_60km_abs_regions(:,:,r) = areas_60km_abs.*(UKregions60 == r);
    areas_60km_frac_regions(:,:,r) = areas_60km_abs_regions(:,:,r)./nansum(nansum(areas_60km_abs_regions(:,:,r)));
end
