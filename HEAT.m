function [] = HEAT(inputs,varargin)
% HEAT v.1.1
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
disp('Running HEAT v.1.1')
disp('-----')

% Set directory paths: essential for running in non-Docker environments
init_HEAT

% Record start time
startt = now;

% For testing purposes, to show the correct data has copied to the Docker
% container (remove this later):
disp(' ')
disp('You are here:')
pwd
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
else
    DAFNI = 0;
end


% Otherwise, if not running on DAFNI:
% Check if inputs is a script, in which case run it
if ischar(inputs)
    if exist(inputs,'file')
        run(inputs)
        
    else
        disp('Error: input script not found');
        return
    end
    
    % Otherwise if inputs is a structure, then adjust the name ready to run HEAT
else
    if ~isstruct(inputs)
        disp('Error: inputs must be a structure');
        return
    end
    if ~isfield(inputs,'ExptName')
        disp('Error: inputs appears to be the wrong structure');
        return
    else
        inputs = inputs;
    end
end


%% Update fields if Environment variables are provided
disp('Updating inputs with environment variables')
% Then overwrite defaults with environment variables if running on DAFNI:
env_BC = getenv('BIASCORR');
env_expn = getenv('EXPNAME');
env_varn = getenv('VARNAME');
env_scen = getenv('SCENARIO');
env_tims = getenv('TIMEPERIOD_S');
env_timl = getenv('TIMEPERIOD_L');

if ~isempty(env_BC)
    env_BC
    if strcmp(env_BC,'y')
        disp('Environment variable found for bias correction option: updating inputs file')
        inputs.BiasCorr = 1;
    end
end

if ~isempty(env_expn)
    disp('Environment variable found for Experiment Name: updating inputs file')
    inputs.ExptName = {env_expn};
%     inputs.ExptName
end
if ~isempty(env_varn)
    disp('Environment variable found for Variable: updating inputs file')
    inputs.Variable = {env_varn};
%     inputs.Variable
end
if ~isempty(env_scen)
    disp('Environment variable found for Scenario: updating inputs file')
    inputs.Scenario = {string(env_scen)};
end

if ~isempty(env_tims)
    disp('Environment variable found for defining time period: updating inputs file')
    inputs.PeriodStart = env_tims;
end

if ~isempty(env_timl)
    inputs.PeriodLength = env_timl;
end
disp(' ')


%% Find which steps of HEAT to run
% Set default to not run steps
runanalysis = 0;
runworkflow = 0;

% Run steps if necessary
if isfield(inputs,'OutputType')
    if strcmp(string(inputs.OutputType),'Analysis')
        runanalysis = 1;
    elseif strcmp(string(inputs.OutputType),'workflow_netCDF')
        runworkflow = 1;
    end
end




%% Set up output directory
% First check the experiment name won't overwrite the DerivedData directory
if strcmp(inputs.ExptName,'DerivedData')
    disp('Cannot call experiment "DerivedData": CANCELLING')
    return
end

% Next, check if output directory already exists
if DAFNI == 0
    Outputdir = [Outputdir,inputs.ExptName,'/'];
end

if exist(Outputdir,'dir')
    
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
    mkdir(Outputdir)
end


% Save input files for future reference
save([Outputdir,'inputs.mat'],'inputs')


%% If an experiment has already been run producing workflow output with this name, delete it
if runworkflow == 1
    
    for v = 1:length(inputs.Variable)
        Variable = char(inputs.Variable(v));
        
        % Check if this file has already been derived:
        froot = Outputdir; % Take the file name...
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

%% Start to load data for analysis

% Load xyz data and regions
load_xyz
load_regions

% Create blank output array if calculating acclimatisation
if isfield(inputs,'MMTpctile')
	reg_acclim = nan(length(inputs.Dataset),13,2); % Dims: models x regions x time slices
    inputs.AnnSummer = 'Annual'; % This only works for annual data
end

%% Go through each UKCP18 simulation
for s = 1:length(inputs.Dataset)
    Dataset = char(inputs.Dataset(s));
    
    if strcmp(inputs.DataType,'UKCP18')
        % Find resolution
        if strcmp(Dataset(1:2),'RC')
            res = '12km/';
            runn = ['run',Dataset(5:6)];
            lats = lat_UK_RCM;
            lons = long_UK_RCM;
            LSM = LSM12;
            UK_area = areas_12km_frac_UK;
            reg_area = areas_12km_frac_regions;
            
        elseif strcmp(Dataset(1:2),'CP')
            res = '2km/';
            runn = ['run',Dataset(5:6)];
            lats = lat_UK_CPM;
            lons = long_UK_CPM;
            LSM = LSM2;
%             reg_area = areas_CPM_frac_regions;
            
        elseif strcmp(Dataset(1:2),'GC')
            res = '60km/';
            runn = ['run',Dataset(5:6)];
            lats = lat_UK_GCM;
            lons = long_UK_GCM;
            LSM = LSM60;
            UK_area = areas_60km_frac_UK;
            reg_area = areas_60km_frac_regions;
            
        elseif strcmp(Dataset(1:2),'CM')
            res = '60km/';
            runn = ['run',Dataset(7:8)];
            lats = lat_UK_GCM;
            lons = long_UK_GCM;
            LSM = LSM60;
            UK_area = areas_60km_frac_UK;
            reg_area = areas_60km_frac_regions;
        end
        
        % Get rid of res sub-directory if using on DAFNI
        if DAFNI == 1
            res = [];
        end
        
        % Update file path and netCDF dims if using bias corr. data
        if inputs.BiasCorr == 1
            if DAFNI == 0
                BCdir = 'BiasCorrected/';
            else 
                BCdir = [];
            end
            % Find location of netCDF data for the required variable
            % Set default domain to load as whole of dataset for T vars
            ncstarts = [1 1 1];
            ncends = [Inf Inf Inf];
            % Dimension of yyyymmdd var for T vars
            datedim = 1;
        else
            BCdir = [];
            % Find location of netCDF data for the required variable
            % Set default domain to load as whole of dataset for T vars
            ncstarts = [1 1 1 1];
            ncends = [Inf Inf Inf Inf];
            % Dimension of yyyymmdd var for T vars
            datedim = 2;
        end
        
        
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
        
        
        % Find list of correct data for this variable
        if strcmp(inputs.Variable,'Tmax')
            % Find what files are available
            if DAFNI == 0
                var = 'tasmax/';
            else 
                var = [];
            end
            % Directory of raw data for each required variable
            vardir = [UKCP18dir,BCdir,res,var,runn,'/'];
            % Find how many files are to be loaded/produced
            files = dir([vardir '*.nc']);
            
        elseif strcmp(inputs.Variable,'Tmean')
            % Find what files are available
            if DAFNI == 0
                var = 'tas/';
                runn = [runn,'/']; % ATKA: Needs repeated for other T vars if this works
            else 
                var = [];
                runn = [];
            end            % Directory of raw data for each required variable
            vardir = [UKCP18dir,BCdir,res,var,runn];
            % Find how many files are to be loaded/produced
            files = dir([vardir '*.nc']);
            
        elseif strcmp(inputs.Variable,'Tmin')
            % Find what files are available
            if DAFNI == 0
                var = 'tasmin/';
            else 
                var = [];
            end            % Directory of raw data for each required variable
            vardir = [UKCP18dir,BCdir,res,var,runn,'/'];
            % Find how many files are to be loaded/produced
            files = dir([vardir '*.nc']);
            
        else % All other heat stress variables are found in  Deriveddir
            var = char(inputs.Variable);
            vardir = Deriveddir;
            % Set default domain to load as whole of dataset for all other vars
            ncstarts = [1 1 1];
            ncends = [Inf Inf Inf];
            % Dimension of yyyymmdd for all other vars
            datedim = 1;
            % Find how many files are to be loaded/produced
            files = dir([vardir,var,'-',inputs.Domain,'-',Dataset,'*.nc']);
        end
        
        % Remove the railing / on the variable name before reading .nc
        var = var(1:end-1);
        
        % Some print outs for testing purposes (remove later):
        disp('This data directory is being accessed:')
        vardir
        disp('These data files are available:')
        ls(vardir)
        disp(' ')
        
        
        %% Go through HadUK-Grid data
        % Set variable, directory and temporal subsetting info for loading
        % HadUK-Grid observations
    elseif strcmp(inputs.DataType,'HadUKGrid')
        % Find resolution
        if strcmp(Dataset(1:2),'12')
            res = '12km/';
            %             runn = ['run',Dataset(5:6)];
            lats = lat_UK_RCM;
            lons = long_UK_RCM;
            LSM = LSM12;
            
        elseif strcmp(Dataset(1:2),'2k')
            res = '2km/';
            %             runn = ['run',Dataset(5:6)];
            lats = lat_UK_CPM;
            lons = long_UK_CPM;
            LSM = LSM2;
            
        elseif strcmp(Dataset(1:2),'60')
            res = '60km/';
            %             runn = ['run',Dataset(5:6)];
            lats = lat_UK_GCM;
            lons = long_UK_GCM;
            LSM = LSM60;
            
        elseif strcmp(Dataset(1:2),'1k')
            res = '1km/';
            %             runn = ['run',Dataset(7:8)];
            lats = lat_UK_HadUK1;
            lons = long_UK_HadUK1;
            LSM = LSM1;
        end
        
        % Find location of netCDF data for the required variable
        % Set default domain to load as whole of dataset for T vars
        ncstarts = [1 1 1];
        ncends = [Inf Inf Inf];
        % Dimension of yyyymmdd var for T vars
        datedim = 2;
        
        if strcmp(inputs.Variable,'Tmax')
            % Find what files are available
            var = 'tasmax';
            % Directory of raw data for each required variable
            vardir = [HadUKdir,var,'/',res];
            % Find how many files are to be loaded/produced
            files = dir([vardir '*.nc']);
            
        elseif strcmp(inputs.Variable,'T')
            % Find what files are available
            var = 'T';
            % Directory of raw data for each required variable
            vardir = Deriveddir;
            files = dir([vardir,'T-',inputs.Domain,'-HadUK-Grid-',Dataset,'*.nc']);
            % Dimension of yyyymmdd for all other vars
            datedim = 1;
            
        elseif strcmp(inputs.Variable,'Tmin')
            % Find what files are available
            var = 'tasmin';
            % Directory of raw data for each required variable
            vardir = [HadUKdir,var,'/',res];
            % Find how many files are to be loaded/produced
            files = dir([vardir '*.nc']);
            
        else % All other heat stress variables are found in  Deriveddir
            var = Variable;
            vardir = Deriveddir;
            % Set default domain to load as whole of dataset for all other vars
            ncstarts = [1 1 1];
            ncends = [Inf Inf Inf];
            % Dimension of yyyymmdd for all other vars
            datedim = 1;
            files = dir([vardir,var,'-',inputs.Domain,'-HadUK-Grid-',Dataset,'*.nc']);
        end
        
        % Find out if temporal subsetting is required (scenarios not
        % available for past observations)
        if isfield(inputs,'PeriodStart')
            TemporalStart = inputs.PeriodStart;
            if isfield(inputs,'PeriodLength')
                TemporalEnd = inputs.PeriodStart + inputs.PeriodLength - 1;
            else
                TemporalEnd = inputs.PeriodStart + 29;
            end
            
            TemporalStart = str2double([num2str(TemporalStart),'0101']);
            TemporalEnd = str2double([num2str(TemporalEnd),'1230']);
            
        end
    end
    
    
    % Assuming data files exist, continue with loading
    if ~isempty(files)
        
        %% Spatially subset data as required
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
        
        
        %% Load only the files required for the temporal subset
        % If specific years are required
        if exist('TemporalStart','var')
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
            file = char([files(i).folder,'/',files(i).name]) % ATKA: Suppress this later if works
            
            % Load temperature for the correct region and concatenate through time if necessary
            if i == startload
                data = double(ncread(file,var,ncstarts,ncends));
                dates = ncread(file,'yyyymmdd');
                times = ncread(file,'time');
                projection_x_coordinate = ncread(file,'projection_x_coordinate',ncstarts(1),ncends(1));
                projection_y_coordinate = ncread(file,'projection_y_coordinate',ncstarts(2),ncends(2));
            else
                data = cat(3,data,double(ncread(file,var,ncstarts,ncends)));
                dates = cat(datedim,dates,ncread(file,'yyyymmdd'));
                times = cat(1,times,ncread(file,'time'));
            end
        end
        
        if datedim == 1
            dates = dates';
        end
        
        
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
                    reg_acclim(s,reg,1) = nansum(nansum(baseline .* UK_area(grid_idx,grid_idy)));
                else
                    reg_acclim(s,reg,1) = nansum(nansum(baseline .* reg_area(grid_idx,grid_idy,reg)));
                end
            end
        end
        
        
        % Pull out the required dates and times
        [data,dates,times] = subset_temporal(data,dates,times,[TemporalStart,TemporalEnd],inputs.AnnSummer);
                
        
        %% Adjust temperature for urban greening
        
        
        
        %% Produce diagnostic data if required (step 2)
        if runanalysis == 1
            
            % Go through each required dataset/simulation/variable
            % Data hierarchy:
            % DataType (e.g. UKCP18, ERA5, CMIP6) -> Dataset (e.g. specific simulation, observation resolution)
            
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
                    reg_acclim(s,reg,2) = nansum(nansum(baseline .* UK_area(grid_idx,grid_idy)));
                else
                    reg_acclim(s,reg,2) = nansum(nansum(baseline .* reg_area(grid_idx,grid_idy,reg)));
                end
            end
        end
        
        
        %% Generate output for other models in workflows
        if runworkflow == 1
            
            % Save as netCDF for HARM to use
            nc_name = [Outputdir,'HEAT-',Dataset,'_',Variable];
            
            xyz.dates = dates;
            xyz.times = times;
            xyz.projection_x_coordinate = projection_x_coordinate;
            xyz.projection_y_coordinate = projection_y_coordinate;
            
            save_HARM_nc(nc_name,data,xyz,'Tmean')
            
%             % csv needs saved for every grid box: go through each lat-long
%             for i = 1:length(data(:,1,1))
%                 for j = 1:length(data(1,:,1))
%                     
%                     % Only save if the point is on land
%                     if LSM(i,j) == 1
%                         
%                         % Set csv output file name
%                         if strcmp(inputs.DataType,'HadUKGrid')
%                             csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_HadUK-Grid-',Dataset,'_',Variable,'.csv'];
%                         elseif strcmp(inputs.DataType,'UKCP18')
%                             csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_UKCP18-',res(1:end-1),'_',Variable,'.csv'];
%                         end
%                         
%                         % Save the data as a .csv file
%                         dlmwrite(csv_name,squeeze(data(i,j,:))','-append','newline','pc','delimiter',',','precision',4);
%                     end
%                 end
%             end
        end
    end
end


%% Save other outputs that might be required later in workflow
% Shift in percentiles for acclimatisation
reg_acclim = reg_acclim(:,:,2) - reg_acclim(:,:,1);
% First, convert to 2D:
reg_acclim_2D = zeros(length(grid_idx),length(grid_idy),12); % lon x lat x model sim
for reg = 1:12
    mask = UKregions12 == reg;
    for sim = 1:12
        reg_acclim_sim = zeros(length(grid_idx),length(grid_idy));
        reg_acclim_sim(mask) = reg_acclim(sim,reg);
        reg_acclim_2D(:,:,sim) = reg_acclim_2D(:,:,sim) + reg_acclim_sim;
        
    end
end
save([Outputdir,'reg_acclim_2D.mat'],'reg_acclim_2D')

% Spatial subset ids
if exist('grid_idx','var')
    save([Outputdir,'grid_idx.mat'],'grid_idx')
end
if exist('grid_idy','var')
    save([Outputdir,'grid_idy.mat'],'grid_idy')
end


%% Finish up
disp(' ')
disp(['HEAT run "',inputs.ExptName,'" complete',])
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))
disp('-----')


