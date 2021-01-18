% generate_region_latlon_area.m
% 
% Like the UK version, this creates the absolute and fractional areas for
% doing spatial averages etc., except for individual regions of England/UK.
% 

%% Load data to start with
load('PreProcessedData/UKregionsERA5.mat')
load('PreProcessedData/UKregions2.mat')
load('PreProcessedData/UKregions12.mat')
load('PreProcessedData/UKregions1.mat')
load('PreProcessedData/UKregions60.mat')
generate_UK_latlon_area

% Regions for reference:
regs = {'Scotland','North East','North West','Yorkshire and the Humber','East Midlands','West Midlands','East of England','Greater London','South East','South West','Wales','Northern Ireland','Isle of Man'};

% Create blank array to fill
areas_ERA5_abs_regions = zeros(length(areas_ERA5_abs(:,1)),length(areas_ERA5_abs(1,:)),12);
areas_ERA5_frac_regions = zeros(length(areas_ERA5_abs(:,1)),length(areas_ERA5_abs(1,:)),12);
areas_RCM_abs_regions = zeros(length(areas_RCM_abs(:,1)),length(areas_RCM_abs(1,:)),12);
areas_RCM_frac_regions = zeros(length(areas_RCM_abs(:,1)),length(areas_RCM_abs(1,:)),12);
areas_GCM_abs_regions = zeros(length(areas_GCM_abs(:,1)),length(areas_GCM_abs(1,:)),12);
areas_GCM_frac_regions = zeros(length(areas_GCM_abs(:,1)),length(areas_GCM_abs(1,:)),12);


for r=1:12
    areas_ERA5_abs_regions(:,:,r) = areas_ERA5_abs.*(UKregionsERA5 == r);
    areas_ERA5_frac_regions(:,:,r) = areas_ERA5_abs_regions(:,:,r)./nansum(nansum(areas_ERA5_abs_regions(:,:,r)));
    areas_RCM_abs_regions(:,:,r) = areas_RCM_abs.*(UKregions12 == r);
    areas_RCM_frac_regions(:,:,r) = areas_RCM_abs_regions(:,:,r)./nansum(nansum(areas_RCM_abs_regions(:,:,r)));
    areas_GCM_abs_regions(:,:,r) = areas_GCM_abs.*(UKregions60 == r);
    areas_GCM_frac_regions(:,:,r) = areas_GCM_abs_regions(:,:,r)./nansum(nansum(areas_GCM_abs_regions(:,:,r)));
end
