function [] = HEAT_step2(inputs,Variable)
% HEAT_step2(inputs,Variable)
%
% Run the second step of HEAT to produce diagnostic outputs for assessing
% heat extremes.
%
% If assessing temperature data, this is loaded from the raw netCDF files.
% If assessing heat stress variables, these are loaded from the derived
% data directory.
%
% Extremes analysis methods that can be assessed are:
%  - Average no. of days exceeding an absolute threshold (inputs.AbsThresh)
%  - Percentile of a metric (inputs.Pctile)
%  - Extreme mean exceeding a given percentile threshold (inputs.ExtremeMeanPctile)
%  - Degree days exceeding an absolute or percentile threshold (inputs.DDa
%   or inputs.DDp respectively)
%
% Output types include a map showing the temporal mean of the extreme
% analysis or a time series showing the spatial mean.

%% Setup
disp(' ')
disp('Running Step 2: producing extremes analysis')
disp('-----')

% Set directory paths
init_HEAT

% Load xyz data
load_xyz

% Set some defaults if required:
% If you haven't specified, you will not get all plots, just MMMs
if ~isfield(inputs,'PlotAll')
    inputs.PlotAll = 0;
end


%% Start analysis
% Set IDs for multi-model mean (if required):
% This will allow different kinds of model simulation to be loaded (e.g.
% GCMs, RCMs) and a MMM taken for each)
MM_id  = zeros(1,length(inputs.Dataset));

% Go through each dataset (e.g. model simulation) requested in inputs
for s = 1:length(inputs.Dataset)
    Dataset = char(inputs.Dataset(s));
    
    %% For UKCP18 data
    for d = 1:length(inputs.DataType)
        DataType = char(inputs.DataType(d));
        if strcmp(DataType,'UKCP18')
            % Define information according to simulation type
            if strcmp(Dataset(1:2),'RC')
                res = '12km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_RCM;
                lons = long_UK_RCM;
                LSM = LSM12;
                areas_abs = areas_12km_abs;
                areas_frac = areas_12km_frac_UK;
                MM_id(s) = 2;
                
            elseif strcmp(Dataset(1:2),'CP')
                res = '2km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_CPM;
                lons = long_UK_CPM;
                LSM = LSM2;
                areas_abs = areas_2km_abs;
                areas_frac = areas_2km_frac_UK;
                MM_id(s) = 3;
                
            elseif strcmp(Dataset(1:2),'GC')
                res = '60km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                areas_abs = areas_60km_abs;
                areas_frac = areas_60km_frac_UK;
                MM_id(s) = 1;
                
            elseif strcmp(Dataset(1:2),'CM')
                res = '60km/';
                runn = ['run',Dataset(7:8)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                areas_abs = areas_60km_abs;
                areas_frac = areas_60km_frac_UK;
                MM_id(s) = 1; % Note: Currently CMIP5 will be included in MMM with UKCP18 GCMs
            end
            
            % Set netCDF information for the required variable:
            % Set default domain to load as whole of dataset for T vars
            ncstarts = [1 1 1 1];
            ncends = [Inf Inf Inf Inf];
            % Dimension of yyyymmdd var for T vars
            datedim = 2;
            
            % Set file locations
            if strcmp(Variable,'Tmax')
                var = 'tasmax';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'T')
                var = 'tas';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'Tmin')
                var = 'tasmin';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,res,var,'/',runn,'/'];
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
            
            
            %% For HadUK-Grid data
        elseif strcmp(DataType,'HadUKGrid')
            % Find resolution
            if strcmp(Dataset(1:2),'12')
                res = '12km/';
                lats = lat_UK_RCM;
                lons = long_UK_RCM;
                LSM = LSM12;
                areas_abs = areas_12km_abs;
                areas_frac = areas_12km_frac_UK;
                
            elseif strcmp(Dataset(1:2),'2k')
                res = '2km/';
                lats = lat_UK_CPM;
                lons = long_UK_CPM;
                LSM = LSM2;
                areas_abs = areas_2km_abs;
                areas_frac = areas_2km_frac_UK;
                
            elseif strcmp(Dataset(1:2),'60')
                res = '60km/';
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                areas_abs = areas_60km_abs;
                areas_frac = areas_60km_frac_UK;
                
            elseif strcmp(Dataset(1:2),'1k')
                res = '1km/';
                lats = lat_UK_HadUK1;
                lons = long_UK_HadUK1;
                LSM = LSM1;
                areas_abs = areas_1km_abs;
                areas_frac = areas_1km_frac_UK;
            end
            
            % Set netCDF information for the required variable:
            % Set default domain to load as whole of dataset for T vars
            ncstarts = [1 1 1];
            ncends = [Inf Inf Inf];
            % Dimension of yyyymmdd var for T vars
            datedim = 2;
            
            % Set file locations
            if strcmp(Variable,'Tmax')
                var = 'tasmax';
                % Directory of raw data for each required variable
                vardir = [HadUKdir,var,'/',res];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'T')
                var = 'T';
                % Directory of raw data for each required variable
                vardir = Deriveddir;
                files = dir([vardir,'T-',inputs.Domain,'-HadUK-Grid-',Dataset,'*.nc']);
                % Dimension of yyyymmdd for all other vars
                datedim = 1;
                
            elseif strcmp(Variable,'Tmin')
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
        
        % Check files exist
        if ~isempty(files)
            
            %% Spatially subset data as required
            % Find corners of requested domain
            if isfield(inputs,'SpatialRange')
                if length(inputs.SpatialRange(1,:)) == 2 % If a lat-long box has been specified to load
                    [lon_id1,lat_id1] = find_location(inputs.SpatialRange(2,1),inputs.SpatialRange(1,1),lons,lats);
                    [lon_id2,lat_id2] = find_location(inputs.SpatialRange(2,2),inputs.SpatialRange(1,2),lons,lats);
                    
                    ncstarts(1) = lon_id1; ncstarts(2) = lat_id1;
                    ncends(1) = 1+lon_id2-lon_id1; ncends(2) = 1+lat_id2-lat_id1;
                    
                    % Create an ID field for subsetting e.g. lat-long, areas etc.
                    grid_idx = lon_id1:lon_id2;
                    grid_idy = lat_id1:lat_id2;
                    
                elseif length(inputs.SpatialRange(1,:)) == 1 % If a specific grid cell has been specified to load
                    [lon_id1,lat_id1] = find_location(inputs.SpatialRange(2,1),inputs.SpatialRange(1,1),lons,lats);
                    
                    ncstarts(1) = lon_id1; ncstarts(2) = lat_id1;
                    ncends(1) = 1; ncends(2) = 1;
                    
                    % Create an ID field for subsetting e.g. lat-long, areas etc.
                    grid_idx = lon_id1;
                    grid_idy = lat_id1;
                    
                end
            else
                grid_idx = 1:length(lats(:,1));
                grid_idy = 1:length(lats(1,:));
            end
            
            
            %% Temporal subsetting
            % Load only the files required for the temporal subset if
            % specific years are required:
            if isfield(inputs,'TemporalRange')
                % Find which files cover the required start and end dates
                for i = 1:length(files)
                    % Find when the netCDF files start and end
                    fstart = files(i).name(end-19:end-12);
                    fend = files(i).name(end-10:end-3);
                    
                    % Find netCDF file that contains required start date
                    if str2double(fstart) <= inputs.TemporalRange(1) && str2double(fend) >= inputs.TemporalRange(1)
                        startload = i;
                    end
                    
                    % Find netCDF file that contains required end date
                    if str2double(fstart) <= inputs.TemporalRange(2) && str2double(fend) >= inputs.TemporalRange(2)
                        endload = i;
                    end
                end
            else
                startload = 1;
                endload = length(files);
                disp('Warning: selected start and/or end date not covered by available data')
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
            
            
            % Some netCDFs have dates saved in a different dimension -
            % convert so they are all the same:
            if datedim == 1
                dates = dates';
            end
            
            
            % Temporally subset to the specific required dates and summer type
            % Pull out the required dates and times
            [data,dates] = subset_temporal(data,dates,inputs.TemporalRange,inputs.AnnSummer);
            
            
            %% Generate map output if required
            for o = 1:length(inputs.OutputType)
                if strcmp(inputs.OutputType(o),'map')
                    
                    % Count number of extreme analysis types requested in the
                    % input file and set an ID for each:
                    n_outputs = 0;
                    EA_type = 0;
                    if isfield(inputs,'Pctile')
                        n_outputs = n_outputs + length(inputs.Pctile);
                        EA_type = cat(2,EA_type,ones(1,length(inputs.Pctile)));
                        EA_1 = 1; % This counter is used later for going through each threshold level (if appropriate)
                    end
                    if isfield(inputs,'ExtremeMeanPctile')
                        n_outputs = n_outputs + 1:length(inputs.ExtremeMeanPctile);
                        EA_type = cat(2,EA_type,ones(1,length(inputs.ExtremeMeanPctile))*2);
                        EA_2 = 1; % This counter is used later for going through each threshold level (if appropriate)
                    end
                    if isfield(inputs,'AbsThresh')
                        n_outputs = n_outputs + length(inputs.AbsThresh);
                        EA_type = cat(2,EA_type,ones(1,length(inputs.AbsThresh))*3);
                        EA_3 = 1; % This counter is used later for going through each threshold level (if appropriate)
                    end
                    
                    % Remove the 0 from the extreme analysis ID array
                    EA_type = EA_type(EA_type>0);
                    
                    % Create output array for all required variables:
                    % 4D: lat x long x model simulation x extreme analysis metrics
                    if s ==1 
                        data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.Dataset),n_outputs);
                    end
                    
                    %% Go through each of the extremes analysis types
                    n_output = 1;
                    for n = 1:n_outputs
                        
                        % Calculate percentile if required
                        if EA_type(n) == 1
                            disp(['Calculating ',num2str(inputs.Pctile(EA_1)),'th percentile over time period'])
                            data_calc(:,:,s,n_output) = prctile(data,inputs.Pctile(EA_1),3);
                            titles2(n_output) = {[Variable, ' ',num2str(inputs.Pctile(EA_1)),'th percentile']};
                            titles1(n_output) = {Dataset};
                            EA_1 = EA_1 + 1; % Update counter so next iteration through loop updates threshold if required
                        end
                        
                        % Calculate extreme mean if required
                        if EA_type(n) == 2
                            disp(['Calculating extreme mean above ',num2str(inputs.ExtremeMeanPctile(EA_2)),'th percentile over time period'])
                            % Calculate percentile threshold
                            Txx_temp = data >= prctile(data,inputs.ExtremeMeanPctile(EA_2),3);
                            Txx = nan(size(Txx_temp));
                            Txx(Txx_temp == 1) = 1;
                            data_Txx = data .* Txx;
                            % Find mean of days exceeding the threshold
                            data_calc(:,:,s,n_output) = squeeze(nanmean(data_Txx,3));
                            titles2(n_output) = {[Variable, ' extreme mean >',num2str(inputs.ExtremeMeanPctile(EA_2)),'th percentile']};
                            titles1(n_output) = {Dataset};
                            EA_2 = EA_2 + 1; % Update counter so next iteration through loop updates threshold if required
                        end
                        
                        % Calculate no. of days exceeding threshold if required
                        if EA_type(n) == 3
                            % Find length of loaded time series (assuming whole summers
                            % have been taken in inputs.TemporalRange)
                            startyr = num2str(inputs.TemporalRange(1));
                            if str2double(startyr(5:6))<=6
                                startyr = str2double(startyr(1:4))-1;
                            else
                                startyr = str2double(startyr(1:4));
                            end
                            
                            endyr = num2str(inputs.TemporalRange(2));
                            if str2double(endyr(5:6))>=7
                                endyr = str2double(endyr(1:4));
                            else
                                endyr = str2double(endyr(1:4))-1;
                            end
                            
                            tslength = endyr-startyr;
                            
                            % Calculate number of days exceeding threshold
                            disp(['Calculating number of days when ',Variable,' exceeds ',num2str(inputs.AbsThresh(EA_3))])
                            data_calc(:,:,s,n_output) = nansum(data>inputs.AbsThresh(EA_3),3)/tslength;
                            titles2(n_output) = {[Variable, '>',num2str(inputs.AbsThresh(EA_3))]};
                            titles1(n_output) = {Dataset};
                            EA_3 = EA_3 + 1; % Update counter so next iteration through loop updates threshold if required
                        end
                        
                        
                        %
                        % Note: Still need to add Degree Day metrics here
                        %
                        
                        
                        %% Plotting maps
                        % If calculating MMM or MMP, only plot all ensemble members if requested
                        if inputs.PlotAll == 1
                            figure
                            UK_subplot(data_calc(:,:,s,n_output) .* LSM(grid_idx,grid_idy),[titles1(n_output),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                            figname = ['map_',num2str(n_output)];
                            saveas(gcf,figname,'png')
                        end
                        
                        
                        % Calculate MMM or MMP
                        % Only do this once all data has been loaded
                        % (i.e. s has reached the end of its loop):
                        if s == length(inputs.Dataset)
                            
                            if isfield(inputs,'MMM')
                                if inputs.MMM == 1
                                    % Only want to mean model simulations, not obs. etc.
                                    % Do this for each model type (MM_ids 1 -> 3):
                                    for id = 1:3
                                        if sum(MM_id == id)>0
                                            data_plot = nanmean(data_calc(:,:,MM_id == id,n_output),3);
                                            
                                            % Plot
                                            figure
                                            title1 = {'MMM'};
                                            UK_subplot(data_plot .* LSM(grid_idx,grid_idy),[title1(1),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                            figname = ['map_MMM_',num2str(id)];
                                            saveas(gcf,figname,'png')
                                        end
                                    end
                                end
                            end
                            
                            if isfield(inputs,'MMP')
                                
                                % Calculate MMP
                                for P = 1:length(inputs.MMP)
                                    % Only want to mean model simulations, not obs. etc.
                                    % Do this for each model type (MM_ids 1 -> 3):
                                    for id = 1:3
                                        if sum(MM_id == id)>0
                                            data_plot = prctile(data_calc(:,:,MM_id == id,n_output),inputs.MMP(P),3);
                                        end
                                        
                                        % Plot
                                        figure
                                        title1 = {['MM ',num2str(inputs.MMP(P)),'th percentile']};
                                        UK_subplot(data_plot .* LSM(grid_idx,grid_idy),[title1(1),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                        figname = ['map_MMP_',num2str(id)];
                                        saveas(gcf,figname,'png')
                                    end
                                end
                            end
                        end % end of MMM/MMP
                        
                    end % end of plotting each extreme analysis type
                    
                    
                %% Generate time series output if required    
                elseif strcmp(inputs.OutputType(o),'timeseries')
                    
                    % Set default averaging length if necessary
                    if ~isfield(inputs,'AveTime')
                        inputs.AveTime = 10;
                    end
                    
                    % Calculate spatial average for selected region to construct time series
                    if isfield(inputs,'Region')
                        % Load region masks
                        load_regions
                        
                        % Create output array for all required variables:
                        % 3D: time x region x model simulation
                        if s == 1
                            data_calc = nan(length(data(1,1,:)),length(inputs.Region),length(inputs.Dataset));
                        end
                        
                        % Go through each requested region
                        for r = 1:length(inputs.Region)
                            for i = 1:12
                                if strcmp(inputs.Region(r),regs(i))
                                    data_calc(:,r,s) = nansum(nansum(data .* areas_reg(:,:,i)));
                                end
                            end
                        end
                    end
                    
                    
                    
                    
                    % Calculate spatial average for selected lat-long box to construct time series
                    if isfield(inputs,'SpatialRange')
                        if length(inputs.SpatialRange(1,:)) == 1
                            
                            % Create output array for all required variables:
                            % 3D: time x region x model simulation
                            if s == 1
                                data_calc = nan(length(data(1,1,:)),length(inputs.Region),length(inputs.Dataset));
                            end
                            
                            % Take time series at selected grid point
                            data_calc(:,r,s) = data(:,:);
                        
                        % Calculate spatial average for selected location to construct time series
                        elseif length(inputs.SpatialRange(1,:)) == 2
                            
                            % Generate area mask for land areas within
                            % selected region:
                            masked_area = LSM(grid_idx,grid_idy) .* areas_abs;
                            masked_area_frac = masked_area ./ nansum(nansum(masked_area));
                            
                            
                            % Create output array for all required variables:
                            % 3D: time x region x model simulation
                            if s == 1
                                data_calc = nan(length(data(1,1,:)),length(inputs.Region),length(inputs.Dataset));
                            end
                            

                            data_calc(:,r,s) = nansum(nansum(data .* masked_area_frac));

                            
                            
                        end
                    end
                    
                    
                end
            end % end of plot style (map vs. time series) selection
        end
    end
end


disp('-----')


