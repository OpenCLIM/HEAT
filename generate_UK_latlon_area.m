% generate_UK_latlon_area.m
%
% Calculates the area of each grid cell in the UK region from the UKCP18,
% ERA5 and HadUK-Grid data for computing spatial means etc.
%


%% Load LSM and grid data
generate_UK_LSM
generate_UK_latlon


%% Area calculation for GCM
if exist('preprocessed/areas_GCM_frac_land.mat','file')
    load preprocessed/areas_GCM_abs.mat
    load preprocessed/areas_GCM_frac.mat
    load preprocessed/areas_GCM_frac_land.mat
    load preprocessed/areas_GCM_frac_UK.mat
    
else
    
    disp('Calculating for GCM resolution')
    
    % Calculate absolute area
    areas_GCM_abs = calc_latlon_area(lat_UK_GCM,long_UK_GCM,'geodarea');
    
    % Calculate fractional area (easier for calculating means later)
    areas_GCM_frac = areas_GCM_abs/nansum(nansum(areas_GCM_abs));
    
    % Calculate fractional area of land areas only
    areas_GCM_frac_land = areas_GCM_abs;
    areas_GCM_frac_land(isnan(LSM60)) = nan;
    areas_GCM_frac_land = areas_GCM_frac_land/nansum(nansum(areas_GCM_abs(~isnan(LSM60))));
    
    % Calculate fractional area of UK land areas only
    areas_GCM_frac_UK = areas_GCM_abs;
    areas_GCM_frac_UK(isnan(datamask60)) = nan;
    areas_GCM_frac_UK = areas_GCM_frac_UK/nansum(nansum(areas_GCM_abs(~isnan(datamask60))));
    
    % Save to speed up in future
    save preprocessed/areas_GCM_abs.mat areas_GCM_abs
    save preprocessed/areas_GCM_frac.mat areas_GCM_frac
    save preprocessed/areas_GCM_frac_land.mat areas_GCM_frac_land
    save preprocessed/areas_GCM_frac_UK.mat areas_GCM_frac_UK
    
end

%% Area calculation for RCM
if exist('preprocessed/areas_RCM_frac_land.mat','file')
    load preprocessed/areas_RCM_abs.mat
    load preprocessed/areas_RCM_frac.mat
    load preprocessed/areas_RCM_frac_land.mat
    load preprocessed/areas_RCM_frac_UK.mat
    
else
    
    disp('Calculating for RCM resolution')

    % Calculate absolute area
    areas_RCM_abs = calc_latlon_area(lat_UK_RCM,long_UK_RCM,'geodarea');
    
    % Calculate fractional area (easier for calculating means later)
    areas_RCM_frac = areas_RCM_abs/nansum(nansum(areas_RCM_abs));
    
    % Calculate fractional area of land areas only
    areas_RCM_frac_land = areas_RCM_abs;
    areas_RCM_frac_land(isnan(LSM12)) = nan;
    areas_RCM_frac_land = areas_RCM_frac_land/nansum(nansum(areas_RCM_abs(~isnan(LSM12))));
    
    % Calculate fractional area of UK land areas only
    areas_RCM_frac_UK = areas_RCM_abs;
    areas_RCM_frac_UK(isnan(datamask12)) = nan;
    areas_RCM_frac_UK = areas_RCM_frac_UK/nansum(nansum(areas_RCM_abs(~isnan(datamask12))));
    
    % Save output
    save preprocessed/areas_RCM_abs.mat areas_RCM_abs
    save preprocessed/areas_RCM_frac.mat areas_RCM_frac
    save preprocessed/areas_RCM_frac_land.mat areas_RCM_frac_land
    save preprocessed/areas_RCM_frac_UK.mat areas_RCM_frac_UK
    
end

%% Area calculation for CPM
if exist('preprocessed/areas_CPM_frac_land.mat','file')
    load preprocessed/areas_CPM_abs.mat
    load preprocessed/areas_CPM_frac.mat
    load preprocessed/areas_CPM_frac_land.mat
    load preprocessed/areas_CPM_frac_UK.mat
    
else
    
    disp('Calculating for CPM resolution')

    % Calculate absolute area
    areas_CPM_abs = calc_latlon_area(lat_UK_CPM,long_UK_CPM,'geodarea');
    
    % Calculate fractional area (easier for calculating means later)
    areas_CPM_frac = areas_CPM_abs/nansum(nansum(areas_CPM_abs));
    
    % Calculate fractional area of land areas only
    areas_CPM_frac_land = areas_CPM_abs;
    areas_CPM_frac_land(isnan(LSM2)) = nan;
    areas_CPM_frac_land = areas_CPM_frac_land/nansum(nansum(areas_CPM_abs(~isnan(LSM2))));
    
    % Calculate fractional area of UK land areas only
    areas_CPM_frac_UK = areas_CPM_abs;
    areas_CPM_frac_UK(isnan(datamask2)) = nan;
    areas_CPM_frac_UK = areas_CPM_frac_UK/nansum(nansum(areas_CPM_abs(~isnan(datamask2))));
    
    save preprocessed/areas_CPM_abs.mat areas_CPM_abs
    save preprocessed/areas_CPM_frac.mat areas_CPM_frac
    save preprocessed/areas_CPM_frac_land.mat areas_CPM_frac_land
    save preprocessed/areas_CPM_frac_UK.mat areas_CPM_frac_UK
    
end

%% Area calculation for ERA5
if exist('preprocessed/areas_ERA5_frac_land.mat','file')
    load preprocessed/areas_ERA5_abs.mat
    load preprocessed/areas_ERA5_frac.mat
    load preprocessed/areas_ERA5_frac_land.mat
    load preprocessed/areas_ERA5_frac_UK.mat
    
else
    
    disp('Calculating for ERA5 resolution')

    % Calculate absolute area
    areas_ERA5_abs = calc_latlon_area(lat_UK_ERA5,long_UK_ERA5,'areaquad');
    
    % Calculate fractional area (easier for calculating means later)
    areas_ERA5_frac = areas_ERA5_abs/nansum(nansum(areas_ERA5_abs));
    
    % Calculate fractional area of land areas only
    areas_ERA5_frac_land = areas_ERA5_abs;
    areas_ERA5_frac_land(isnan(LSMERA5)) = nan;
    areas_ERA5_frac_land = areas_ERA5_frac_land/nansum(nansum(areas_ERA5_abs(~isnan(LSMERA5))));
    
    % Calculate fractional area of UK land areas only
    areas_ERA5_frac_UK = areas_ERA5_abs;
    areas_ERA5_frac_UK(isnan(datamaskERA5)) = nan;
    areas_ERA5_frac_UK = areas_ERA5_frac_UK/nansum(nansum(areas_ERA5_abs(~isnan(datamaskERA5))));
    
    % Save output
    save preprocessed/areas_ERA5_abs.mat areas_ERA5_abs
    save preprocessed/areas_ERA5_frac.mat areas_ERA5_frac
    save preprocessed/areas_ERA5_frac_land.mat areas_ERA5_frac_land
    save preprocessed/areas_ERA5_frac_UK.mat areas_ERA5_frac_UK
    
end

%% Area calculation for HadUK-Grid 1km
if exist('preprocessed/areas_HadUK1_frac_UK.mat','file')
    load preprocessed/areas_HadUK1_abs.mat
    load preprocessed/areas_HadUK1_frac.mat
    load preprocessed/areas_HadUK1_frac_UK.mat
    
else
    
    disp('Calculating for HadUK-Grid 1km resolution')

    % Calculate absolute area
    areas_HadUK1_abs = calc_latlon_area(lat_UK_HadUK1,long_UK_HadUK1,'geodarea');
    
    % Calculate fractional area (easier for calculating means later)
    areas_HadUK1_frac = areas_HadUK1_abs/nansum(nansum(areas_HadUK1_abs));
    
    % Calculate fractional area of UK land areas only
    areas_HadUK1_frac_UK = areas_HadUK1_abs;
    areas_HadUK1_frac_UK(isnan(datamask1)) = nan;
    areas_HadUK1_frac_UK = areas_HadUK1_frac_UK/nansum(nansum(areas_HadUK1_abs(~isnan(datamask1))));
    
    % Save output
    save preprocessed/areas_HadUK1_abs.mat areas_HadUK1_abs
    save preprocessed/areas_HadUK1_frac.mat areas_HadUK1_frac
    save preprocessed/areas_HadUK1_frac_UK.mat areas_HadUK1_frac_UK
    
end

%% Area calculation for HadUK-Grid 12km
if exist('preprocessed/areas_HadUK12_frac_UK.mat','file')
    
    load preprocessed/areas_HadUK12_abs.mat
    load preprocessed/areas_HadUK12_frac.mat
    load preprocessed/areas_HadUK12_frac_UK.mat
    
else
    
    disp('Calculating for HadUK-Grid 12km resolution')
    
    % Calculate absolute area
    areas_HadUK12_abs = calc_latlon_area(lat_UK_HadUK12,long_UK_HadUK12,'geodarea');
    
    
    % Calculate fractional area (easier for calculating means later)
    areas_HadUK12_frac = areas_HadUK12_abs/nansum(nansum(areas_HadUK12_abs));
    
    % Calculate fractional area of UK land areas only
    areas_HadUK12_frac_UK = areas_HadUK12_abs;
    areas_HadUK12_frac_UK = HadUK2UKCP18(areas_HadUK12_frac_UK); % Correct grid
    areas_HadUK12_frac_UK(isnan(datamask12)) = nan;
    areas_HadUK12_frac_UK = areas_HadUK12_frac_UK/nansum(nansum(areas_HadUK12_abs(~isnan(datamask12))));
    
    % Save output
    save preprocessed/areas_HadUK12_abs.mat areas_HadUK12_abs
    save preprocessed/areas_HadUK12_frac.mat areas_HadUK12_frac
    save preprocessed/areas_HadUK12_frac_UK.mat areas_HadUK12_frac_UK
    
end

%% Area calculation for HadUK-Grid 60km
if exist('preprocessed/areas_HadUK60_frac_UK.mat','file')
    
    load preprocessed/areas_HadUK60_abs.mat
    load preprocessed/areas_HadUK60_frac.mat
    load preprocessed/areas_HadUK60_frac_UK.mat
    
else
    
    disp('Calculating for HadUK-Grid 60km resolution')
    
    % Calculate absolute area
    areas_HadUK60_abs = calc_latlon_area(lat_UK_HadUK60,long_UK_HadUK60,'geodarea');
    
    % Calculate fractional area (easier for calculating means later)
    areas_HadUK60_frac = areas_HadUK60_abs/nansum(nansum(areas_HadUK60_abs));
    
    % Calculate fractional area of land areas only
    areas_HadUK60_frac_UK = areas_HadUK60_abs;
    areas_HadUK60_frac_UK(isnan(datamask60)) = nan;
    areas_HadUK60_frac_UK = areas_HadUK60_frac_UK/nansum(nansum(areas_HadUK60_abs(~isnan(datamask60))));
    
    % Save output
    save preprocessed/areas_HadUK60_abs.mat areas_HadUK60_abs
    save preprocessed/areas_HadUK60_frac.mat areas_HadUK60_frac
    save preprocessed/areas_HadUK60_frac_UK.mat areas_HadUK60_frac_UK
    
end

