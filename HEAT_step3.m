function [] = HEAT_step3(inputs,DataType,Variable)
% HEAT_step3(inputs,DataType,Variable)
%
% Run the third step of HEAT to output workflow data. This is data that is
% used by other OpenCLIM models in workflows on DAFNI. Currently, the only
% workflow that HEAT feeds into is the ARCADIA heat-mortality model. Other
% models can be added as required.
%
% It uses the same loading process as HEAT_step2, whereby each dataset
% (e.g. simulation) is loaded in turn either from raw or derived data for
% the correct spatial and temporal domain. There is an option to output
% regional means for the ARCADIA model ? if a region is specified then
% whole spatial domain is loaded before taking the regional mean.
%
% ARCADIA currently only requires daily mean temperature or sWBGT as far
% as I know. Daily max. and daily min. temperatures can be exported also.
%
% TO DO: Check if ARCADIA always requires a full annual cycle of data. Katie's
% document on Teams suggests so, in which case the AnnSummer subsetting
% should be removed and the calendar for UKCP18 models needs converted from
% 360 to 365 days (I can use the ISIMIP2b method for this).

%% Generate ARCADIA-ready csv if requested
disp(' ')
disp('Running Step 3: producing workflow output for ARCADIA')
disp('-----')

% Find out if temporal subsetting is required
if isfield(inputs,'TemporalRange')
    TemporalStart = inputs.TemporalRange(1);
    TemporalEnd = inputs.TemporalRange(2);
end


if strcmp(inputs.WorkflowOutput,'ARCADIA')
    
    % Set directory paths
    init_HEAT
    
    % Load xyz data
    load_xyz
    
    %% Go through each simulation
    for s = 1:length(inputs.Dataset)
        Dataset = char(inputs.Dataset(s));
        
        %% Go through UKCP18 data
        if strcmp(DataType,'UKCP18')
            % Find resolution
            if strcmp(Dataset(1:2),'RC')
                res = '12km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_RCM;
                lons = long_UK_RCM;
                LSM = LSM12;
                
            elseif strcmp(Dataset(1:2),'CP')
                res = '2km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_CPM;
                lons = long_UK_CPM;
                LSM = LSM2;
                
            elseif strcmp(Dataset(1:2),'GC')
                res = '60km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                
            elseif strcmp(Dataset(1:2),'CM')
                res = '60km/';
                runn = ['run',Dataset(7:8)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
            end
            
            
            % Update file path and netCDF dims if using bias corr. data
            if inputs.BiasCorr == 1
                BCdir = 'BiasCorrected/';
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
            
            
            % Find if scenario is required, and if so, set temporal start and
            % end:
            if isfield(inputs,'Scenario')
                % Load the years each scenario reaches a warming level
                load('PreProcessedData/tas_GCM_glob_thresh_arr.mat')
                % Look up the correct model run from the table
                modelslist = {'run01','run04','run05','run06','run07','run08','run09','run10','run11','run12','run13','run15'};
                modelid = find(contains(modelslist,runn));
                
                % Read the start year of period
                if strcmp(inputs.Scenario,'past')
                    TemporalStart = 1990;
                elseif strcmp(inputs.Scenario,'1.5')
                    TemporalStart = tas_GCM_glob_thresh_arr(1,modelid);
                elseif strcmp(inputs.Scenario,'2.0')
                    TemporalStart = tas_GCM_glob_thresh_arr(2,modelid);
                elseif strcmp(inputs.Scenario,'3.0')
                    TemporalStart = tas_GCM_glob_thresh_arr(4,modelid);
                elseif strcmp(inputs.Scenario,'4.0')
                    TemporalStart = tas_GCM_glob_thresh_arr(6,modelid);
                end
                
                % Get end year
                TemporalEnd = TemporalStart+29;
                
                % Convert years to correct format
                TemporalStart = str2double([num2str(TemporalStart),'0101']);
                TemporalEnd = str2double([num2str(TemporalEnd),'1230']);
            end

            
            
            if strcmp(Variable,'Tmax')
                % Find what files are available
                var = 'tasmax';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,BCdir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'Tmean')
                % Find what files are available
                var = 'tas';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,BCdir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'Tmin')
                % Find what files are available
                var = 'tasmin';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,BCdir,res,var,'/',runn,'/'];
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
                % Find how many files are to be loaded/produced
                files = dir([vardir,var,'-',inputs.Domain,'-',Dataset,'*.nc']);
            end
            
            vardir
            ls(vardir)
            
            
            %% Go through HadUK-Grid data
        elseif strcmp(DataType,'HadUKGrid')
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
            
            if strcmp(Variable,'Tmax')
                % Find what files are available
                var = 'tasmax';
                % Directory of raw data for each required variable
                vardir = [HadUKdir,var,'/',res];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'T')
                % Find what files are available
                var = 'T';
                % Directory of raw data for each required variable
                vardir = Deriveddir;
                files = dir([vardir,'T-',inputs.Domain,'-HadUK-Grid-',Dataset,'*.nc']);
                % Dimension of yyyymmdd for all other vars
                datedim = 1;
                
            elseif strcmp(Variable,'Tmin')
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
        end
        
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
                file = [files(i).folder,'/',files(i).name];
                
                % Load temperature for the correct region and concatenate through time if necessary
                if i == startload
                    data = double(ncread(file,var,ncstarts,ncends));
                    dates = ncread(file,'yyyymmdd');
                else
                    data = cat(3,data,double(ncread(file,var,ncstarts,ncends)));
                    dates = cat(datedim,dates,ncread(file,'yyyymmdd'));
                end
            end
            
            if datedim == 1
                dates = dates';
            end
            
            
            %% Temporally subset to the specific required dates and summer type
            if ~isfield(inputs,'AnnSummer')
                inputs.AnnSummer = 'Annual';
            end
            
            % Pull out the required dates and times
            [data,dates] = subset_temporal(data,dates,[TemporalStart,TemporalEnd],inputs.AnnSummer);
            
            
            %% Output ARCADIA data for grid cell or regional mean
            % If interested in every single grid cell
            if ~isfield(inputs,'OutputRegion')
                
                
                % csv needs saved for every grid box: go through each lat-long
                for i = 1:length(data(:,1,1))
                    for j = 1:length(data(1,:,1))
                        
                        % Only save if the point is on land
                        if LSM(i,j) == 1
                            
                            % Set csv output file name
                            if strcmp(DataType,'HadUKGrid')
                                csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_HadUK-Grid-',Dataset,'_',Variable,'.csv'];
                            elseif strcmp(DataType,'UKCP18')
                                csv_name = [Outputdir,'/',inputs.ExptName,'/',num2str(i),'_',num2str(j),'_UKCP18-',res(1:end-1),'_',Variable,'.csv'];
                            end
                            
                            
                            %                         % If csv has been created already
                            %                         if exist(csv_name,'file')
                            %                             % Load existing csv
                            %                             data_csv = csvread(csv_name);
                            %                             % Extract individual grid cell for current dataset
                            %                             data_csv_temp = squeeze(data(i,j,:))';
                            %                             % Concatenate this with the existing csv and resave
                            %                             data_csv = cat(1,data_csv,data_csv_temp);
                            %                             csvwrite(csv_name,data_csv)
                            %                         else
                            %                             % Otherwise extract individual grid cell for current
                            %                             % dataset and create new csv
                            %                             data_csv = squeeze(data(i,j,:))';
                            %                             csvwrite(csv_name,data_csv)
                            %                         end
                            
                            % Save the data as a .csv file
                            dlmwrite(csv_name,squeeze(data(i,j,:))','-append','newline','pc','delimiter',',','precision',4);
                        end
                    end
                end
                
                
            else % Otherwise take a regional mean
                regmean = calc_reg_mean(data,inputs.Region);
                
                
                % Set csv output file name
                csv_name = [Outputdir,'/',inputs.ExptName,'/',inputs.Region,'_',inputs.ExptName,'_',Variable,'.csv'];
                
%                 % If csv has been created already
%                 if exist(csv_name,'file')
%                     % Load existing csv
%                     data_csv = csvread(csv_name);
%                     
%                     % Concatenate this with the existing csv and resave
%                     data_csv = cat(1,data_csv,regmean');
%                     csvwrite(csv_name,data_csv)
%                 else
%                     % Otherwise create new csv
%                     csvwrite(csv_name,regmean')
%                 end
                
                % Save the data as a .csv file
                dlmwrite(csv_name,regmean','-append','newline','pc','delimiter',',','precision',4);
                
            end
        end
        
    end
end
end

