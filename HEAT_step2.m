function [] = HEAT_step2(inputs,DataType,Variable)
% HEAT_step2(inputs,Dataset,Variable)
%
% Run the second step of HEAT to produce diagnostic outputs for assessing
% heat extremes.
%
% If assessing temperature data, this is loaded from the raw netCDF files.
% If assessing heat stress variables, these are loaded from the derived
% data directory.
%

%% Setup
disp(' ')
disp('Running Step 2: producing extremes analysis')
disp('-----')

% Set directory paths
init_HEAT

% Load xyz data
load_xyz

% Set some defaults if required:
% If you haven't specified, you will get all plots, not just MMMs
if ~isfield(inputs,'PlotAll')
    inputs.PlotAll = 1;
end

%% Go through each simulation
for s = 1:length(inputs.Dataset)
    Dataset = char(inputs.Dataset(s));
    
    %% Go through UKCP18 data
    for d = 1:length(inputs.DataType)
        DataType = char(inputs.DataType(d));
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
            
            % Find location of netCDF data for the required variable
            % Set default domain to load as whole of dataset for T vars
            ncstarts = [1 1 1 1];
            ncends = [Inf Inf Inf Inf];
            % Dimension of yyyymmdd var for T vars
            datedim = 2;
            
            if strcmp(Variable,'Tmax')
                % Find what files are available
                var = 'tasmax';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'T')
                % Find what files are available
                var = 'tas';
                % Directory of raw data for each required variable
                vardir = [UKCP18dir,res,var,'/',runn,'/'];
                % Find how many files are to be loaded/produced
                files = dir([vardir '*.nc']);
                
            elseif strcmp(Variable,'Tmin')
                % Find what files are available
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
            if isfield(inputs,'TemporalRange')
                % Find which files cover the required start and end dates
                for i = 1:length(files)
                    % Find when the netCDF files start and end
                    fstart = files(i).name(end-19:end-12);
                    fend = files(i).name(end-10:end-3);
                    
                    % Find netCDF file that contains required start
                    if str2double(fstart) <= inputs.TemporalRange(1) && str2double(fend) >= inputs.TemporalRange(1)
                        startload = i;
                    end
                    
                    % Find netCDF file that contains required end
                    if str2double(fstart) <= inputs.TemporalRange(2) && str2double(fend) >= inputs.TemporalRange(2)
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
            % Pull out the required dates and times
            [data,dates] = subset_temporal(data,dates,inputs.TemporalRange,inputs.AnnSummer);
            
            
            %% Calculate percentile if required
            if isfield(inputs,'Pctile')
                
                % Create empty array for output
                data_pctiles = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.Pctile));
                
                % Calculate given percentile
                for p = 1:length(inputs.Pctile)
                    disp(['Calculating ',num2str(inputs.Pctile(p)),'th percentile over time period'])
                    data_pctiles(:,:,p) = prctile(data,inputs.Pctile(p),3);
                    
                    % Plot as required
                    for o = 1:length(inputs.OutputType)
                        OutputType = string(inputs.OutputType(o));
                        
                        % If calculating MMM or MMP, only plot all ensemble members if requested
                        if inputs.PlotAll ~= 0
                            
                            % Generate map if requested
                            if strcmp(OutputType,'map')
                                figure
                                UK_subplot(data_pctiles(:,:,p) .* LSM(grid_idx,grid_idy),[Dataset,' ',Variable, ' ',num2str(inputs.Pctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                
                            end
                        end
                        
                    end
                    
                end
                
                % Store data for MMM/MMP if required and calculate
                if inputs.MMM == 1 || ~isempty(inputs.MMP)
                    % Only want to mean model simulations, not obs. etc.
                    if strcmp(Dataset(1),'1') || strcmp(Dataset(1),'2') || strcmp(Dataset(1),'6')
                        data_pctiles_mod = nan(size(data_pctiles));
                    else
                        data_pctiles_mod = data_pctiles;
                    end
                    
                    
                    % If this is the first simulation loaded, make new array
                    if s == 1
                        data_all = data_pctiles_mod;
                        
                        % Otherwise, add the current simulation to the array
                    elseif s > 1 && s < length(inputs.Dataset)
                        data_all = cat(ndims(data_pctiles)+1,data_all,data_pctiles_mod);
                        
                        % And if this is the final simulation, calculate the MMM or MMP
                    elseif s == length(inputs.Dataset)
                        data_all = cat(ndims(data_pctiles)+1,data_all,data_pctiles_mod);
                        
                        % Calculate MMM or MMP
                        if inputs.MMM == 1
                            data_plot = nanmean(data_all,ndims(data_pctiles)+1);
                        elseif ~isempty(inputs.MMP)
                            data_plot = prctile(data_all,inputs.MMP,ndims(data_pctiles)+1);
                        end
                        
                        % Go through each precentile to plot
                        for p = 1:length(inputs.Pctile)
                            % Plot as required
                            for o = 1:length(inputs.OutputType)
                                OutputType = string(inputs.OutputType(o));
                                
                                % Generate map if requested
                                if strcmp(OutputType,'map')
                                    figure
                                    UK_subplot(data_plot(:,:,p) .* LSM(grid_idx,grid_idy),['MMM ',Variable, ' ',num2str(inputs.Pctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                    
                                end
                                
                            end
                        end
                        
                    end
                    
                    
                end
                
                
                
            end
            
            %% Calculate extreme mean if required
            if isfield(inputs,'ExtremeMeanPctile')
                
                % Create empty array for output
                data_ExMean = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.ExtremeMeanPctile));
                
                for p = 1:length(inputs.ExtremeMeanPctile)
                    disp(['Calculating extreme mean above ',num2str(inputs.ExtremeMeanPctile(p)),'th percentile over time period'])
                    
                    % Calculate percentile threshold
                    Txx_temp = data >= prctile(data,inputs.ExtremeMeanPctile(p),3);
                    Txx = nan(size(Txx_temp));
                    Txx(Txx_temp == 1) = 1;
                    data_Txx = data .* Txx;
                    % Find mean of days exceeding the threshold
                    data_ExMean(:,:,p) = squeeze(nanmean(data_Txx,3));
                    
                    % Plot as required
                    for o = 1:length(inputs.OutputType)
                        OutputType = string(inputs.OutputType(o));
                        
                        % If calculating MMM or MMP, only plot all ensemble members if requested
                        if inputs.PlotAll ~= 0
                            
                            % Generate map if requested
                            if strcmp(OutputType,'map')
                                figure
                                UK_subplot(data_ExMean(:,:,p) .* LSM(grid_idx,grid_idy),[Dataset,' ',Variable, ' extreme mean >',num2str(inputs.ExtremeMeanPctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                
                            end
                        end
                        
                    end
                    
                end
                
                
                % Store data for MMM/MMP if required and calculate
                if inputs.MMM == 1 || ~isempty(inputs.MMP)
                    % Only want to mean model simulations, not obs. etc.
                    if strcmp(Dataset(1),'1') || strcmp(Dataset(1),'2') || strcmp(Dataset(1),'6')
                        data_ExMean_mod = nan(size(data_ExMean));
                    else
                        data_ExMean_mod = data_ExMean;
                    end
                    
                    
                    % If this is the first simulation loaded, make new array
                    if s == 1
                        data_all = data_ExMean_mod;
                        
                        % Otherwise, add the current simulation to the array
                    elseif s > 1 && s < length(inputs.Dataset)
                        data_all = cat(ndims(data_ExMean)+1,data_all,data_ExMean_mod);
                        
                        % And if this is the final simulation, calculate the MMM or MMP
                    elseif s == length(inputs.Dataset)
                        data_all = cat(ndims(data_ExMean)+1,data_all,data_ExMean_mod);
                        
                        % Calculate MMM or MMP
                        if inputs.MMM == 1
                            data_plot = nanmean(data_all,ndims(data_ExMean)+1);
                        elseif ~isempty(inputs.MMP)
                            data_plot = prctile(data_all,inputs.MMP,ndims(data_ExMean)+1);
                        end
                        
                        % Go through each precentile to plot
                        for p = 1:length(inputs.ExtremeMeanPctile)
                            % Plot as required
                            for o = 1:length(inputs.OutputType)
                                OutputType = string(inputs.OutputType(o));
                                
                                % Generate map if requested
                                if strcmp(OutputType,'map')
                                    figure
                                    UK_subplot(data_plot(:,:,p) .* LSM(grid_idx,grid_idy),['MMM ',Variable, ' extreme mean >',num2str(inputs.ExtremeMeanPctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                                    
                                end
                                
                            end
                        end
                        
                    end
                    
                    
                end
            end
        end
        
        %         %% Store data for MMM if required and calculate
        %         if inputs.MMM == 1 || ~isempty(inputs.MMP)
        %             % If this is the first simulation loaded, make new array
        %             if s == 1
        %                 data_all = data;
        %
        %             % If this is the final simulation, calculate the MMM or MMP
        %             elseif s == length(inputs.Dataset)
        %                 data_all = cat(ndims(data)+1,data_all,data);
        %
        %                 % Calculate MMM or MMP
        %                 if inputs.MMM == 1
        %                     data = nanmean(data_all,ndims(data)+1);
        %
        %                 elseif ~isempty(inputs.MMP)
        %                     data = prctile(data_all,inputs.MMP,ndims(data)+1);
        %
        %                 end
        %
        %                 % Otherwise, add the current simulation to the array
        %             else
        %                 data_all = cat(ndims(data)+1,data_all,data);
        %             end
        %         end
        
        
        
        
        
        
        %     %% Calculate and produce analysis outputs
        %     % Calculate percentile if required
        %     if isfield(inputs,'Pctile')
        %
        %         % Create empty array for output
        %         data_pctiles = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.Pctile));
        %
        %         for p = 1:length(inputs.Pctile)
        %             disp(['Calculating ',num2str(inputs.Pctile(p)),'th percentile over time period'])
        %
        %             data_pctiles(:,:,p) = prctile(data,inputs.Pctile(p),3);
        %
        %             % Plot as required
        %             for o = 1:length(inputs.OutputType)
        %                 OutputType = string(inputs.OutputType(o));
        %
        %                 % Generate map if requested
        %                 if strcmp(OutputType,'map')
        %                     figure
        %                     UK_subplot(data_pctiles(:,:,p),[Dataset,' ',Variable, ' ',num2str(inputs.Pctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy))
        %                     colorbar()
        %                     % Save as required
        %
        %                 end
        %
        %             end
        %
        %
        %
        %         end
        %     end
        %
        %     % Calculate extreme mean if required
        %     if isfield(inputs,'ExtremeMeanPctile')
        %
        %         % Create empty array for output
        %         data_ExMean = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.ExtremeMeanPctile));
        %
        %         for p = 1:length(inputs.ExtremeMeanPctile)
        %             disp(['Calculating extreme mean above ',num2str(inputs.ExtremeMeanPctile(p)),'th percentile over time period'])
        %
        %
        %             Txx_temp = data >= prctile(data,inputs.ExtremeMeanPctile(p),3);
        %             Txx = nan(size(Txx_temp));
        %             Txx(Txx_temp == 1) = 1;
        %             data_Txx = data .* Txx;
        %             data_ExMean(:,:,p) = squeeze(nanmean(data_Txx,3));
        %
        %             % Plot as required
        %             for o = 1:length(inputs.OutputType)
        %                 OutputType = string(inputs.OutputType(o));
        %
        %                 % Generate map if requested
        %                 if strcmp(OutputType,'map')
        %                     figure
        %                     UK_subplot(data_ExMean(:,:,p),[Dataset,' ',Variable, ' extreme mean >',num2str(inputs.Pctile(p)),'th percentile'],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy))
        %                     colorbar()
        %                     % Save as required
        %
        %                 end
        %
        %             end
        %
        %             % Save as required
        %
        %         end
        %     end
        
        
        
        %     %% Calculate the MMM or MMP if required
        %     if inputs.MMM == 1
        %         data = nanmean(data_all,ndims(data)+1);
        %
        %     elseif ~isempty(inputs.MMP)
        %         data = prctile(data_all,inputs.MMP,ndims(data)+1);
        %
        %     end
    end
end


% % Generate requested output
% disp('Producing output')
% for o = 1:length(inputs.OutputType)
%     OutputType = string(inputs.OutputType(o));
%
%     %     % Generate map if requested
%     %     if strcmp(OutputType,'map')
%     %         figure
%     %         UK_subplot(data_ExMean,[Dataset,' ',Variable, ' extreme mean >95th percentile'])
%     %         %                             figure
%     %         %                             UK_subplot(data_pctiles,[Dataset,' ',Variable, ' 95th percentile'])
%     %     end
%
%
%
%
% end
disp('-----')


