function [] = HEAT(inputs,Climatedata,Urbandirin)
% HEAT v.2.0
%
% Run the OpenCLIM Heat Extremes Analysis Toolbox (HEAT). This function
% does an inital setup and check of file directories, then goes through
% each required dataset, loading, processing and generating output
% accordingly.
%
% Input files, in the form of structures, can be created using the format
% provided in launch_file_TEMPLATE.m or interactively using
% generate_launch_file.m. Alternatively the inputs can be defined directly
% from a script (e.g. launch_file_TEMPLATE.m) if running HEAT app
% standalone.
%
% inputs = information regarding data which is to be loaded, including
% variables, dataset to use, temporal and spatial domain, experiment name
% and whether to save derived variables.
%
% Coded by A.T. Kennedy-Asser, University of Bristol, 2022.
% Contact: alan.kennedy@bristol.ac.uk
%

%% Initialise
disp('Running HEAT v.2.0')
disp('-----')

% Set directory paths: essential for running in non-Docker environments
init_HEAT

% Record start time
startt = now;

% For testing purposes, to show the correct data has copied to the Docker
% container (remove this later):
if strcmp(pwd,'/code')
    disp('Running in Docker container with these files:')
    ls
    disp('-----')
    disp(' ')

end

%% Check inputs
% If running HEAT as standalone app, inputs should be included in a
% seperate inputs script, which can be adapted from input_files_TEMPLATE.m.
% If running HEAT through the MATLAB GUI, then it is possible to have run
% input_files_TEMPLATE.m first, in which case the structures this script
% creates can be input.

% If running on DAFNI, no 'inputs' will have been passed directly to HEAT.m
% Load the default DAFNI template:
if ~exist('inputs','var')
    disp('No inputs file passed to Docker: running default for DAFNI')
    input_files_DAFNI
    disp(' ')
    
    % Create output directory
    mkdir(Climatedirout)
    
% Otherwise, running locally for testing:
else
    % Find out if ClimateData has been passed as an argument and, if so, is it
    % a file or a directory of files
    if exist(Climatedata,'file')
        Climatedirin = Climatedata;
    else
        disp('Valid climate data missing: CANCELLING')
        disp('-----')
        return
    end
    
    % Setup output directory if testing locally:
    Climatedirout = [Climatedirout,inputs.ExptName,'/'];
    
    if exist(Climatedirout,'dir')
        
        % Check if overwriting has been authorised in input file
        if ~isfield(inputs,'OverwriteExpt')
            disp('Permission to overwrite unclear in input file: CANCELLING')
            disp('-----')
            % Cancel the run
            return
        else
            
            % If it does, check whether or not to overwrite
            disp('Warning: existing experiment exists with this name')
            if inputs.OverwriteExpt == 1
                disp('Overwriting enabled in input file: Re-running experiment')
                disp('-----')
                
            else
                disp('Overwriting disabled in input file: CANCELLING')
                disp('-----')
                % If not overwriting, cancel the run
                return
            end
            
        end
        
    else
        % Create output directory
        mkdir(Climatedirout)
    end
end

% Save input files for future reference
save([Climatedirout,'inputs.mat'],'inputs')


%% Find which steps of HEAT to run
% Set default to not run steps
runexmean = 0;
runDD = 0;
runanalysis = 0;
runworkflow = 0;
runabsext = 0;
runperext = 0;


% Run steps if necessary
if isfield(inputs,'OutputType')    
    if strcmp(string(inputs.OutputType),'ExtremeMean') 
        runexmean = 1;
    elseif strcmp(string(inputs.OutputType),'DD66')
        runDD = 1;
    elseif strcmp(string(inputs.OutputType),'AbsoluteExtremes')
        runabsext = 1;
    elseif strcmp(string(inputs.OutputType),'PercentileExtremes')
        runperext = 1;
    elseif strcmp(string(inputs.OutputType),'Heatwave exposure (HE)') % NEEDS ADDED
        runanalysis = 1;
    elseif strcmp(string(inputs.OutputType),'NetCDF')
        runworkflow = 1;
    end
end


%% If an experiment has already been run producing workflow output with this name, delete it
if runworkflow == 1
    for v = 1:length(inputs.Variable)
        Variable = char(inputs.Variable(v));
        
        % Check if this file has already been derived:
        froot = Climatedirout; % Take the file name...
        files = dir([froot '*.nc']); % Then check if any files exist with this root
        
        % If so, delete
        if ~isempty(files)
            disp('HARM-ready .nc files already exist: deleting these.')
            for f = 1:length(files)
                file = [files(f).folder,'/',files(f).name];
                delete(file)
            end
        end
    end
end


%% Load other default data and setup blank output arrays
% Load xyz data and regions
load_xyz
load_regions

% Create blank output array if calculating acclimatisation
if isfield(inputs,'MMTpctile')
    reg_acclim = nan(13,2); % Dims: regions x time slices
    inputs.AnnSummer = 'Annual'; % This only works for annual data
end


%% Adjust temperature for urban greening
if exist('Urbandirin','var')
    % Find list of files to load
    urbfiles = dir([Urbandirin '*.asc']);
    
    % Assuming data files exist, continue with loading
    if isempty(urbfiles)
        disp('No urban development data to load: no change made to UHI intensity')
        disp('-----')
        disp(' ')
    else
        disp('The following urban development data is available to be loaded:')
        ls([Urbandirin '*.asc'])
        disp('-----')
        if exist([Urbandirin 'out_cell_dev-12km-sum.asc'],'file')
            [baseline_urb,RefMat]= arcgridread([Urbandirin,'out_cell_dev-12km-sum.asc']);
        else
            [baseline_urb,RefMat]= arcgridread([Urbandirin,'out_cell_dev.asc']);
        end
        
        disp('Urban data loaded')
        
        % Find existing development
        dev_old = baseline_urb == 1;
        % Find new development
        dev_new = baseline_urb >= 1;
        save([Climatedirout,'dev_old.mat'],'dev_old')
        save([Climatedirout,'dev_new.mat'],'dev_new')
        
        % If on the correct grid for the RCM, don't change res
        if length(baseline_urb(:,1)) == 82 && length(baseline_urb(1,:)) == 112
            
            % Find change in urban isation
            urb_change = dev_new - dev_old;
            save([Climatedirout,'urb_change.mat'],'urb_change')
            dev_all = dev_new;
            
        % Regrid only if necessary    
        else
            % Prepare a lat/lon grid
            [nrows,ncols,~]=size(baseline_urb);
            [row,col]=ndgrid(1:nrows,1:ncols);
            [lat_urb,lon_urb]=pix2latlon(RefMat,row,col);
            [lats_urb,lons_urb] = os2ll(lon_urb,lat_urb);
            disp('Urban data regridded to lat-long')
            
            % Aggregate to 2km
            regrid_res = 20; % FOR TESTING
            
            x = floor(length(baseline_urb(:,1))/regrid_res);
            y = floor(length(baseline_urb(1,:))/regrid_res);
            dev_old_1km = zeros(x,y);
            dev_new_1km = zeros(x,y);
            lat_urb_1km = nan(x,y);
            lon_urb_1km = nan(x,y);
            for i = 1:x
                for j = 1:y
                    ii = 1+(i-1)*regrid_res:regrid_res+(i-1)*regrid_res;
                    jj = 1+(j-1)*regrid_res:regrid_res+(j-1)*regrid_res;
                    
                    dev_old_1km(i,j) = nanmean(nanmean(dev_old(ii,jj)));
                    dev_new_1km(i,j) = nanmean(nanmean(dev_new(ii,jj)));
                    lat_urb_1km(i,j) = nanmean(nanmean(lats_urb(ii,jj)));
                    lon_urb_1km(i,j) = nanmean(nanmean(lons_urb(ii,jj)));
                end
            end
            disp('Urban data aggregated to 2km')
            
            % Find change in urban area
            urb_change = dev_new_1km - dev_old_1km;
            save([Climatedirout,'urb_change.mat'],'urb_change')
            
            % Re-grid to RCM grid
            disp('Starting re-gridding')
            dev_new_interp = griddata(lon_urb_1km,lat_urb_1km,dev_new_1km,long_UK_RCM,lat_UK_RCM,'linear');
            urb_change_interp = griddata(lon_urb_1km,lat_urb_1km,urb_change,long_UK_RCM,lat_UK_RCM,'linear');
            disp('Re-gridding complete')
            
            urb_change_interp_raw = urb_change_interp;
            save([Climatedirout,'urb_change_interp_raw.mat'],'urb_change_interp_raw')
            urb_change_interp(urb_change_interp_raw > 1) = 1;
            save([Climatedirout,'urb_change_interp.mat'],'urb_change_interp')
            
            % Give standard name for use later
            urb_change = urb_change_interp;
            dev_all = dev_new_interp;
        end
        
        % Adjust temperature based on increased UHI intensity
        UHI_I = str2double(inputs.UHI_I); % Value based upon offline analysis (load_urban_fraction.m), default = 2, plausible range ~ 1.5 - 3. Possibly include option to change this in future.
        
        UHI_adjustment = urb_change * UHI_I;
        UHI_adjustment(isnan(UHI_adjustment)) = 0;
        %         data = data + UHI_adjustment;
        disp('UHI adjjustment complete')
        save([Climatedirout,'UHI_adjustment.mat'],'UHI_adjustment')
        csvwrite([Climatedirout,'UHI_adjustment.csv'],UHI_adjustment)
        disp('-----')
        disp(' ')

    end
end


%% Urban greening parameterisation
if isfield(inputs,'Greening')
    if exist('dev_all','var')
        GreenEffect = dev_all *double(inputs.Greening);
        GreenEffect(isnan(GreenEffect)) = 0;
    end
end


%% Find list of files to load and define parameters for input climate dataset
files = dir([Climatedirin '*.nc']);

% Assuming data files exist, continue with loading
if isempty(files)
    disp('No valid climate data to load: CANCELLING')
    return
else
    disp('The following climate data netCDFs are being loaded:')
    ls([Climatedirin '*.nc'])
    disp('-----')    
    disp(' ')
end

% Find if using one of the UKCP18 RCM ensemble members
if isfield(inputs,'Dataset')
    Dataset = char(inputs.Dataset(1));
    if strcmp(Dataset(1:2),'RC')
        runn = ['run',Dataset(5:6)];
    end
else
    Dataset = 'data';
end


% Load one file as an example to check if data is in 3D
file = ([files(1).folder,'/',files(1).name]);
datatest = double(ncread(file,char(inputs.Variable)));
s1 = size(datatest);
ncstarts = [1 1 1];
ncends = [Inf Inf Inf];

if ndims(datatest) > 3
    disp('Input netCDF has too many dimensions (more than 3): CANCELLING')
    return
end

% Also find out orientation of date strings (assuming there is a long timeseries)
datestest = ncread(file,'yyyymmdd');
s2 = size(datestest);
if s2(1) > s2(2)
    datedim = 1;
else
    datedim = 2;
end
    
    % Find resolution
    disp('Checking model resolution:')
if s1(1) == 82 && s1(2) == 112
    disp('Input data is on UKCP18 RCM grid')
    lats = lat_UK_RCM;
    lons = long_UK_RCM;
    LSM = LSM12;
    UK_area = areas_12km_frac_UK;
    reg_area = areas_12km_frac_regions;
    
    subsetting = 1;
    averaging = 1;
    
elseif s1(1) == 17 && s1(2) == 23
    disp('Input data is on UKCP18 GCM grid')
    lats = lat_UK_GCM;
    lons = long_UK_GCM;
    LSM = LSM60;
    UK_area = areas_60km_frac_UK;
    reg_area = areas_60km_frac_regions;
    
    subsetting = 1;
    averaging = 1;
    
else
    disp('Input data is on unknown grid: Area subsetting and regional percentile acclimatisation disabled.')
    LSM = ones(s1(1),s1(2));
    
    subsetting = 0;
    averaging = 0;
end
disp('-----')
disp(' ')



% % Update file path and netCDF dims if using bias corr. data
% if inputs.BiasCorr == 1
%     % Find location of netCDF data for the required variable
%     % Set default domain to load as whole of dataset for T vars
%     ncstarts = [1 1 1];
%     ncends = [Inf Inf Inf];
%     % Dimension of yyyymmdd var for T vars
%     datedim = 1;
% else
%     % Find location of netCDF data for the required variable
%     % Set default domain to load as whole of dataset for T vars
%     ncstarts = [1 1 1 1];
%     ncends = [Inf Inf Inf Inf];
%     % Dimension of yyyymmdd var for T vars
%     datedim = 2;
% end


%% Set up spatial subset if required
% Find corners of requested domain
if subsetting == 1
    if isfield(inputs,'SpatialRange')
        if length(inputs.SpatialRange(1,:)) == 2 % lat-long box specified to load
            [lon_id1,lat_id1] = find_location(inputs.SpatialRange(2,1),inputs.SpatialRange(1,1),lons,lats);
            [lon_id2,lat_id2] = find_location(inputs.SpatialRange(2,2),inputs.SpatialRange(1,2),lons,lats);
            
            ncstarts(1) = lon_id1; ncstarts(2) = lat_id1;
            ncends(1) = 1+lon_id2-lon_id1; ncends(2) = 1+lat_id2-lat_id1;
            
            % Create an ID field for subsetting e.g. lat-long, areas etc.
            grid_idx = lon_id1:lon_id2;
            grid_idy = lat_id1:lat_id2;
            
        elseif length(inputs.SpatialRange(1,:)) == 1 % specific grid cell specified
            [lon_id1,lat_id1] = find_location(inputs.SpatialRange(2,1),inputs.SpatialRange(1,1),lons,lats);
            
            ncstarts(1) = lon_id1; ncstarts(2) = lat_id1;
            ncends(1) = 1; ncends(2) = 1;
            
            % Create an ID field for subsetting e.g. lat-long, areas etc.
            grid_idx = lon_id1;
            grid_idy = lat_id1;
            
        end
        
    elseif isfield(inputs,'Region')
        % TO DO: Insert code for selecting specific regions
    end
end


%% Setup temporal subsetting
% Find out if temporal subsetting is required or if scenario is
% required, and if so, set temporal start and end:
if isfield(inputs,'PeriodStart')
    TemporalStart = inputs.PeriodStart;
    if isfield(inputs,'PeriodLength')
        TemporalEnd = inputs.PeriodStart + inputs.PeriodLength - 1;
    else
        TemporalEnd = inputs.PeriodStart + 29;
    end
    
    TemporalStart = str2double([num2str(TemporalStart),'0101']);
    TemporalEnd = str2double([num2str(TemporalEnd),'1230']);
    
elseif isfield(inputs,'Scenario')
    % Load the years each scenario reaches a warming level
    load('PreProcessedData/tas_GCM_glob_thresh_arr_arnell.mat')
    % TO DO: add option so users can upload their own time slice
    % info
    
    % Look up the correct model run from the table
    modelslist = {'run01','run04','run05','run06','run07','run08','run09','run10','run11','run12','run13','run15'};
    modelid = find(contains(modelslist,runn));
    
    % Read the start year of period
    if strcmp(inputs.Scenario,'past')
        TemporalStart = 1990;
    elseif strcmp(inputs.Scenario,'s1.5')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(1,modelid);
    elseif strcmp(inputs.Scenario,'s2.0')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(2,modelid);
    elseif strcmp(inputs.Scenario,'s3.0')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(4,modelid);
    elseif strcmp(inputs.Scenario,'s4.0')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(6,modelid);
    end
    
    % Set max time period to 2050-2079 (end of RCM simulations)
    if TemporalStart > 2050
        TemporalStart = 2050;
    end
    
    % Get end year
    TemporalEnd = TemporalStart+29;
    
    % Convert years to correct format
    TemporalStart = str2double([num2str(TemporalStart),'0101']);
    TemporalEnd = str2double([num2str(TemporalEnd),'1230']);
else
    disp('No time period defined: STOPPING')
    return
end


%% Load only the files required for the temporal/spatial subset
% If specific years are required
if exist('PeriodStart','var')
    % Find which files cover the required start and end dates
    for i = 1:length(files)
        % Find when the netCDF files start and end
        fstart = files(i).name(end-19:end-12);
        fend = files(i).name(end-10:end-3);
        
        % Find netCDF file that contains required start
        if str2double(fstart) <= TemporalStart && str2double(fend) >= TemporalStart
            startload = i;
        end
        
        % Find netCDF file that contains required end
        if str2double(fstart) <= TemporalEnd && str2double(fend) >= TemporalEnd
            endload = i;
        end
    end
else
    startload = 1;
    endload = length(files);
end

% Load all of the files between the start and end file
for i = startload:endload
    
    % File name
    file = ([files(i).folder,'/',files(i).name]);
    
    % Load temperature for the correct region and concatenate through time if necessary
    if i == startload
        data = double(ncread(file,char(inputs.Variable),ncstarts,ncends));
        dates = ncread(file,'yyyymmdd');
        times = ncread(file,'time');
        projection_x_coordinate = ncread(file,'projection_x_coordinate',ncstarts(1),ncends(1));
        projection_y_coordinate = ncread(file,'projection_y_coordinate',ncstarts(2),ncends(2));
    else
        data = cat(3,data,double(ncread(file,char(inputs.Variable),ncstarts,ncends)));
        dates = cat(datedim,dates,ncread(file,'yyyymmdd'));
        times = cat(1,times,ncread(file,'time'));
    end
end

if datedim == 1
    dates = dates';
end


%% Tidy up temporally
% Tidy up the calendar dates:

% Update TemporalEnd if there is data for 31st December
if strcmp(string(dates(5:8,end)'),'1231')
    TemporalEnd = TemporalEnd + 1;
end

% Remove leap year days
dates_no_ly = string(dates(5:8,:)');
keep_dates = ~strcmp(dates_no_ly,'0229');
times = times(keep_dates);
dates = dates(:,keep_dates);
data = data(:,:,keep_dates);


% Temporally subset to the specific required dates and summer type
if ~isfield(inputs,'AnnSummer')
    inputs.AnnSummer = 'Annual';
end

if ~exist('grid_idx','var')
    grid_idx = [1:length(projection_x_coordinate)];
end
if ~exist('grid_idy','var')
    grid_idy = [1:length(projection_y_coordinate)];
end

% If required, calculate acclimatisation baseline
if averaging == 1
    if isfield(inputs,'MMTpctile')
        baseline = prctile(data(:,:,1:7300),inputs.MMTpctile,3); % This takes 1981-2000 as the baseline
        for reg = 1:13
            if reg == 13
                reg_acclim(reg,1) = nansum(nansum(baseline .* UK_area(grid_idx,grid_idy)));
            else
                reg_acclim(reg,1) = nansum(nansum(baseline .* reg_area(grid_idx,grid_idy,reg)));
            end
        end
    end
end

% Pull out the required dates and times
[data,dates,times] = subset_temporal(data,dates,times,[TemporalStart,TemporalEnd],inputs.AnnSummer);


% Adjust for UHI, if available
if exist('UHI_adjustment','var')
    
    datarange = [nanmin(nanmin((mean(data,3).*LSM))) nanmax(nanmax((mean(data,3).*LSM)))];
    
    % Produce some output to sanity check
    UKave = nansum(nansum(mean(data,3) .* UK_area(grid_idx,grid_idy)));
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(mean(data,3).*LSM,['Data, no UHI (mean = ',num2str(UKave),')'],Climatedirout,lat_UK_RCM,long_UK_RCM,datarange)
    end
    
    data = data + UHI_adjustment;
    UKave = nansum(nansum(mean(data,3) .* UK_area(grid_idx,grid_idy)));
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(mean(data,3).*LSM,['Data + UHI (mean = ',num2str(UKave),')'],Climatedirout,lat_UK_RCM,long_UK_RCM,datarange)
    end
end

% Adjust for Greening effect, if available
if exist('GreenEffect','var')
    data = data + GreenEffect;
    UKave = nansum(nansum(mean(data,3) .* UK_area(grid_idx,grid_idy)));
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(mean(data,3).*LSM,['Data + UHI (mean = ',num2str(UKave),')'],Climatedirout,lat_UK_RCM,long_UK_RCM,datarange)
    end
end


%% Calculate extreme mean if required
if runexmean == 1
    % Set default if necessary
    if ~isfield(inputs,'ExtremeMeanPctile')
        inputs.ExtremeMeanPctile = 95; 
        disp('Calculating extreme mean for default 95th percentile')
        disp('For reproduction of Kennedy-Asser et al. 2022, use Tmax summer only data')
        disp('-----')
    end
    
    % Find xth percentile for each grid point
    exmeanthresh = prctile(data,inputs.ExtremeMeanPctile,3);
    
    % Find days that do not exceed xth percentile
    nonexdays = data < exmeanthresh;
    
    % Copy data then remove non-extreme days
    exdays = data;
    exdays(nonexdays) = nan;
    
    % Calculate extreme mean
    exmean = nanmean(exdays,3);
    
    % Save output
%     save([Climatedirout,'exmean.mat'],'exmean')
    dlmwrite([Climatedirout,'exmean.csv'],exmean, 'delimiter', ',', 'precision', '%i')
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(exmean.*LSM,'Extreme mean',Climatedirout,lat_UK_RCM,long_UK_RCM)
    end

end


%% Calculate DD66 degree day metric if required
if runDD == 1
    % Set default if necessary
    if ~isfield(inputs,'DD')
        inputs.DD = 66; 
        disp('Calculating Degree Days for default 66th percentile')
        disp('For reproduction of Kennedy-Asser et al. 2022, use Tmean summer only data')
        disp('-----')
    end
    
    % Find xth percentile for each grid point
    DDthresh = prctile(data,inputs.DD,3);
    
    % Find days that do not exceed xth percentile
    nonexdays = data < DDthresh;
    
    % Copy data then remove non-extreme days
    exdays = data - DDthresh;
    exdays(nonexdays) = nan;
    
    % If it hasn't been defined already, assume the period length is 30
    % years
    if ~isfield(inputs,'PeriodLength')
        inputs.PeriodLength = 30;
    end
    
    % Calculate extreme mean
    DDx = nansum(exdays,3) / inputs.PeriodLength;
    
    % Save output
%     save([Climatedirout,'DDx.mat'],'DDx')
    dlmwrite([Climatedirout,'DDx.csv'],DDx, 'delimiter', ',', 'precision', '%i')
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(DDx.*LSM,'Degree Days',Climatedirout,lat_UK_RCM,long_UK_RCM)
    end
    
end


%% Calculate number of days exceeding absolute extreme value
if runabsext == 1
    disp('Days above absolute threshold')
    isfield(inputs,'AbsThresh')
    % Set default if necessary
    if ~isfield(inputs,'AbsThresh')
        inputs.AbsThresh = 25; 
        disp('Calculating number of days exceeding 25 degC (default)')
        disp('-----')
    end
    
    % If it hasn't been defined already, assume the period length is 30
    % years
    if ~isfield(inputs,'PeriodLength')
        inputs.PeriodLength = 30;
    end
    
    % Calculate average number of days > x �C per year
    AbsExt = nansum(data > inputs.AbsThresh,3) / inputs.PeriodLength;
    
    % Save output
    dlmwrite([Climatedirout,'AbsExt.csv'],AbsExt, 'delimiter', ',', 'precision', '%i')
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(AbsExt.*LSM,['Number of days exceeding ',num2str(inputs.AbsThresh),' degC'],Climatedirout,lat_UK_RCM,long_UK_RCM)
    end
    
end


%% Calculate number of days exceeding percentile extreme value
if runperext == 1
    disp('Days above absolute threshold')
    isfield(inputs,'PercentileThresh')
    % Set default if necessary
    if ~isfield(inputs,'PercentileThresh')
        inputs.PercentileThresh = 95; 
        disp('Calculating number of days exceeding 95th percentile (default)')
        disp('-----')
    end
    
    
    Thresh = prctile(data,inputs.PercentileThresh,3);
    
    % If it hasn't been defined already, assume the period length is 30
    % years
    if ~isfield(inputs,'PeriodLength')
        inputs.PeriodLength = 30;
    end
    
    % Calculate average number of days > xth percentile per year
    PerExt = nansum(data > Thresh,3) / inputs.PeriodLength;
    
    % Save output
    dlmwrite([Climatedirout,'Threshold.csv'],Thresh, 'delimiter', ',', 'precision', '%i')
    dlmwrite([Climatedirout,'NDays.csv'],PerExt, 'delimiter', ',', 'precision', '%i')
    if ~isfield(inputs,'SpatialRange')
    figure
    UK_subplot(PerExt.*LSM,['Number of days exceeding ',num2str(inputs.PercentileThresh),'th percentile'],Climatedirout,lat_UK_RCM,long_UK_RCM)
    end
    
end



%% Calculate acclimatisation as shift in regional percentile
if averaging == 1
    if isfield(inputs,'MMTpctile')
        baseline = prctile(data,inputs.MMTpctile,3);
        for reg = 1:13
            if reg == 13
                reg_acclim(reg,2) = nansum(nansum(baseline .* UK_area(grid_idx,grid_idy)));
            else
                reg_acclim(reg,2) = nansum(nansum(baseline .* reg_area(grid_idx,grid_idy,reg)));
            end
        end
    end
end


%% Generate output for other models in workflows
if runworkflow == 1
    
    % Save as netCDF for HARM to use
    nc_name = [Climatedirout,'HEAT-',inputs.ExptName,'-',Dataset,'_',Variable];
    
    xyz.dates = dates;
    xyz.times = times;
    xyz.projection_x_coordinate = projection_x_coordinate;
    xyz.projection_y_coordinate = projection_y_coordinate;
    
    % Check dim sizes
    disp('Data dimension sizes = ')
    size(data)
    disp('X dimension sizes = ')
    size(xyz.projection_x_coordinate)
    disp('Y dimension sizes = ')
    size(xyz.projection_y_coordinate)
    disp('Time dimension sizes = ')
    length(data(1,1,:))
    
    save_HARM_nc(nc_name,data,xyz,char(inputs.Variable))
end



%% Save other outputs that might be required later in workflow
if isfield(inputs,'MMTpctile')
    % Shift in percentiles for acclimatisation
    reg_acclim = reg_acclim(:,2) - reg_acclim(:,1);
    % First, convert to 2D:
    reg_acclim_2D = zeros(length(grid_idx),length(grid_idy)); % lon x lat
    for reg = 1:12
        mask = UKregions12 == reg;
        reg_acclim_sim = zeros(length(grid_idx),length(grid_idy));
        reg_acclim_sim(mask) = reg_acclim(reg);
        reg_acclim_2D(:,:) = reg_acclim_2D(:,:) + reg_acclim_sim;
        
    end
    
    save([Climatedirout,'reg_acclim_2D.mat'],'reg_acclim_2D')
end

% Spatial subset ids
if exist('grid_idx','var')
    save([Climatedirout,'grid_idx.mat'],'grid_idx')
end
if exist('grid_idy','var')
    save([Climatedirout,'grid_idy.mat'],'grid_idy')
end


%% Finish up
disp(' ')
disp(['HEAT run "',inputs.ExptName,'" complete',])
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))
disp('-----')
close all
