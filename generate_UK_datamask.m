% generate_UK_datamask.m
% 
% Script to load in some HadUK-Grid data to create a datamask (i.e. showing
% areas where data is available). This is for all land areas of the UK so
% is used to inform the generation of the LSMs in generate_UK_LSM.m.

%% 'Standard' datamasks (i.e. on native HadUK grids)
% Load if already calculated
if exist('preprocessed/datamask1.mat','file')
    load('preprocessed/datamask1.mat')
    load('preprocessed/datamask12.mat')
    load('preprocessed/datamask60.mat')
    
% Otherwise generate    
else
    % 1km resolution
    [datamask1] = load_HadUK_data('tasmax','1km','JJA',1);
    datamask1 = double(~isnan(datamask1(:,:,1)));
    
    % 12km resolution
    [datamask12] = load_HadUK_data('tasmax','12km','JJA',1);
    datamask12 = ~isnan(datamask12(:,:,1));
    datamask12 = HadUK2UKCP18(datamask12); % Correct the grid
    
    % 60km resolution
    [datamask60] = load_HadUK_data('tasmax','60km','JJA',1);
    datamask60 = double(~isnan(datamask60(:,:,1)));
    
    % Make sure values are NaN or 1
    datamask1(datamask1 == 0) = nan;
    datamask12(datamask12 == 0) = nan;
    datamask60(datamask60 == 0) = nan;

    % Save
    save('preprocessed/datamask1.mat','datamask1')
    save('preprocessed/datamask12.mat','datamask12')
    save('preprocessed/datamask60.mat','datamask60')

end

%% Interpolated datamasks for other resolutions
% Load if already calculated
if exist('preprocessed/datamask2.mat','file')
    load('preprocessed/datamask2.mat')
    load('preprocessed/datamaskERA5.mat')
    
% Otherwise generate    
else
    
    % Load lat-long data for interpolating to other resolutions
    generate_UK_latlon

    % Interpolate 1km to 2.2km
    disp('Interpolating to 2.2 km grid')
    datamask2 = griddata(long_UK_HadUK1,lat_UK_HadUK1,double(datamask1),long_UK_CPM,lat_UK_CPM);
    datamask2(datamask2<0.5) = nan;
    datamask2(datamask2>=0.5) = 1;  
    
    % Interpolate 12km to 0.5° resolution
    disp('Interpolating to ERA5 grid')
    datamaskERA5 = griddata(long_UK_HadUK1,lat_UK_HadUK1,double(datamask1),long_UK_ERA5,lat_UK_ERA5);
    datamaskERA5(datamaskERA5<0.5) = nan;
    datamaskERA5(datamaskERA5>=0.5) = 1;

    % Make sure values are NaN or 1
    datamask2(datamask2 == 0) = nan;
    datamaskERA5(datamaskERA5 == 0) = nan;
    
    % Save
    save('preprocessed/datamask2.mat','datamask2')
    save('preprocessed/datamaskERA5.mat','datamaskERA5')
    
end

%% Remove Shetland land points on edge of domain in datamask2
datamask2a = datamask2;
datamask2a(:,606) = nan;

%% Create smaller ERA5 datamask for UK only
datamaskERA5a = datamaskERA5(49:88,42:83);
