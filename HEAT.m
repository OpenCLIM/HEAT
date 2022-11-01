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
    disp(' ')
    disp('-----')
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
    DAFNI = 1; % A helpful flag for later on
    input_files_DAFNI
    disp(' ')
    
    % Create output directory
    mkdir(Climatedirout)
    
    % Otherwise, running locally for testing:
else
    DAFNI = 0;
    
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

% Run steps if necessary
if isfield(inputs,'OutputType')
    if strcmp(string(inputs.OutputType),'Extreme mean') % NEEDS ADDED
        runexmean = 1;
    elseif strcmp(string(inputs.OutputType),'DD66') % NEEDS ADDED
        runDD = 1;
    elseif strcmp(string(inputs.OutputType),'Absolute extremes') % NEEDS ADDED
        runanalysis = 1;
    elseif strcmp(string(inputs.OutputType),'Percentile extremes') % NEEDS ADDED
        runanalysis = 1;
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


%% Define parameters for input climate dataset
Dataset = char(inputs.Dataset(1));

% Find resolution
if strcmp(Dataset(1:2),'RC')
    runn = ['run',Dataset(5:6)];
    lats = lat_UK_RCM;
    lons = long_UK_RCM;
    LSM = LSM12;
    UK_area = areas_12km_frac_UK;
    reg_area = areas_12km_frac_regions;
end


% Update file path and netCDF dims if using bias corr. data
if inputs.BiasCorr == 1
    % Find location of netCDF data for the required variable
    % Set default domain to load as whole of dataset for T vars
    ncstarts = [1 1 1];
    ncends = [Inf Inf Inf];
    % Dimension of yyyymmdd var for T vars
    datedim = 1;
else
    % Find location of netCDF data for the required variable
    % Set default domain to load as whole of dataset for T vars
    ncstarts = [1 1 1 1];
    ncends = [Inf Inf Inf Inf];
    % Dimension of yyyymmdd var for T vars
    datedim = 2;
end


%% Set up spatial subset if required
% Find corners of requested domain
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
    elseif strcmp(inputs.Scenario,'1.5')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(1,modelid);
    elseif strcmp(inputs.Scenario,'2.0')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(2,modelid);
    elseif strcmp(inputs.Scenario,'3.0')
        TemporalStart = tas_GCM_glob_thresh_arr_arnell(4,modelid);
    elseif strcmp(inputs.Scenario,'4.0')
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
end


%% Find list of files to load
files = dir([Climatedirin '*.nc']);

% Assuming data files exist, continue with loading
if isempty(files)
    disp('No valid climate data to load: CANCELLING')
    return
else
    disp('The following climate data netCDFs are being loaded:')
    ls([Climatedirin '*.nc'])
    disp('-----')
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

% Pull out the required dates and times
[data,dates,times] = subset_temporal(data,dates,times,[TemporalStart,TemporalEnd],inputs.AnnSummer);


%% Adjust temperature for urban greening
if exist('Urbandirin','var')
    % Find list of files to load
    urbfiles = dir([Urbandirin '*.asc']);
    
    % Assuming data files exist, continue with loading
    if isempty(urbfiles)
        disp('No urban development data to load: no change made to UHI intensity')
    else
        disp('The following urban development data is available to be loaded:')
        ls([Urbandirin '*.asc'])
        disp('-----')
        [baseline_urb,RefMat]= arcgridread([Urbandirin,'out_cell_dev.asc']);
        
        % Prepare a lat/lon grid
        [nrows,ncols,~]=size(baseline_urb);
        [row,col]=ndgrid(1:nrows,1:ncols);
        [lat_urb,lon_urb]=pix2latlon(RefMat,row,col);
        [lats_urb,lons_urb] = os2ll(lon_urb,lat_urb);
        
        % Find existing development
        dev_old = baseline_urb == 1;
        % Find new development
        dev_new = baseline_urb >= 1;
        
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
        
        % Find change in urban area
        urb_change = dev_new_1km - dev_old_1km;
        
        % Re-grid to RCM grid
        dev_old_interp = griddata(lon_urb_1km,lat_urb_1km,dev_old_1km,long_UK_RCM,lat_UK_RCM,'linear');
        urb_change_interp = griddata(lon_urb_1km,lat_urb_1km,urb_change,long_UK_RCM,lat_UK_RCM,'linear');

        % Adjust temperature based on increased UHI intensity
        UHI_I = 2; % Value based upon offline analysis (load_urban_fraction.m), plausible range ~ 1.5 - 3. Possibly include option to change this in future.
        
        UHI_adjustment = urb_change_interp * UHI_I;
        UHI_adjustment(isnan(UHI_adjustment)) = 0;
        data = data + UHI_adjustment;
    end
end


%% Calculate extreme mean if required
if runexmean == 1
    % Set default if necessary
    if ~isfield(inputs,'ExtremeMeanPctile')
        inputs.ExtremeMeanPctile = 95; 
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
    save('exmean.mat','exmean')
end


%% Calculate DD66 degree day metric if required
if runDD66 == 1
    % Set default if necessary
    if ~isfield(inputs,'DD')
        inputs.DD = 66; % Set default if necessary
    end
    
    
end


%% Produce diagnostic data if required (step 2)
if runanalysis == 1
    
    % Load each required variable
    for v = 1:length(inputs.Variable)
        Variable = char(inputs.Variable(v));
        
        % Run Step 2: Extremes analysis
        if runanalysis == 1
            HEAT_step2(inputs,Variable)
        end
    end
end


%         %% Extract only heatwave days if required
%         % Extract days depending on threshold definition:
%         % If an absolute threshold map has been provided
%         if exist('Thresholds','dir')
%             disp('Extracting heatwave days based on absolute threshold map')
%             load('Thresholds/threshmap.mat')
%             % Note: this currently assumes the threshold map will be
%             % provided as a raster .mat file. Flexibility to provide this
%             % in other ways will need to be added later.
%             [HWdays,numdays,avelength,numevents] = extract_consecHWdays(data,threshmap-10,inputs.Duration,'absolute');
%             data(HWdays==0) = nan;
%             disp(' ')
%
%         % If using a percentile threshold across all of UK
%         elseif isfield(inputs,'PctThresh')
%             disp('Extracting heatwave days based on percentile threshold')
%             [HWdays,numdays,avelength,numevents] = extract_consecHWdays(data,inputs.PctThresh,inputs.Duration,'percentile');
%             disp(' ')
%
%         % If using a single absolute threshold across all of UK
%         elseif isfield(inputs,'AbsThresh')
%             disp('Extracting heatwave days based on absolute threshold')
%             [HWdays,numdays,avelength,numevents] = extract_consecHWdays(data,inputs.AbsThresh,inputs.Duration,'absolute');
%             disp(' ')
%
%         end


%% Calculate acclimatisation as shift in regional percentile
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


%% Generate output for other models in workflows
if runworkflow == 1
    
    % Save as netCDF for HARM to use
    nc_name = [Climatedirout,'HEAT-',Dataset,'_',Variable];
    
    xyz.dates = dates;
    xyz.times = times;
    xyz.projection_x_coordinate = projection_x_coordinate;
    xyz.projection_y_coordinate = projection_y_coordinate;
    
    save_HARM_nc(nc_name,data,xyz,'tas')
end



%% Save other outputs that might be required later in workflow
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


