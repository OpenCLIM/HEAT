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
% If you haven't specified, you will get all plots, not just MMMs
if ~isfield(inputs,'PlotAll')
    inputs.PlotAll = 1;
end


%% Go through each simulation
% Set IDs for MMM (if required)
MM_id  = zeros(1,length(inputs.Dataset));

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
                MM_id(s) = 2;
                
            elseif strcmp(Dataset(1:2),'CP')
                res = '2km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_CPM;
                lons = long_UK_CPM;
                LSM = LSM2;
                MM_id(s) = 3;
                
            elseif strcmp(Dataset(1:2),'GC')
                res = '60km/';
                runn = ['run',Dataset(5:6)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                MM_id(s) = 1;
                
            elseif strcmp(Dataset(1:2),'CM')
                res = '60km/';
                runn = ['run',Dataset(7:8)];
                lats = lat_UK_GCM;
                lons = long_UK_GCM;
                LSM = LSM60;
                MM_id(s) = 1;
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
            
            
            %% Generate map output if required
            if strcmp(inputs.OutputType,'map')
                
                % Count number of extreme analysis types
                n_outputs = 0;
                EA_type = 0;
                if isfield(inputs,'Pctile')
                    n_outputs = n_outputs + length(inputs.Pctile);
                    EA_type = cat(2,EA_type,ones(1,length(inputs.Pctile)));
                    EA_1 = 1;
                end
                if isfield(inputs,'ExtremeMeanPctile')
                    n_outputs = n_outputs + 1:length(inputs.ExtremeMeanPctile);
                    EA_type = cat(2,EA_type,ones(1,length(inputs.ExtremeMeanPctile))*2);
                    EA_2 = 1;
                end
                if isfield(inputs,'AbsThresh')
                    n_outputs = n_outputs + length(inputs.AbsThresh);
                    EA_type = cat(2,EA_type,ones(1,length(inputs.AbsThresh))*3);
                    EA_3 = 1;
                end
                
                EA_type = EA_type(EA_type>0);
                
                % Create output array for all required variables:
                % 4D: lat x long x model simulation x extremes metric
                data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.Dataset),n_outputs);
                
                % Go through each of the extremes analysis types
                n_output = 1;
                for n = 1:n_outputs
                    
                    %% Calculate percentile if required
                    if EA_type(n) == 1
                        disp(['Calculating ',num2str(inputs.Pctile(EA_1)),'th percentile over time period'])
                        data_calc(:,:,s,n_output) = prctile(data,inputs.Pctile(EA_1),3);
                        titles2(n_output) = {[Variable, ' ',num2str(inputs.Pctile(EA_1)),'th percentile']};
                        titles1(n_output) = {Dataset};
                        EA_1 = EA_1 + 1;
                    end
                    
                    %% Calculate extreme mean if required
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
                        EA_2 = EA_2 + 1;
                    end
                    
                    %% Calculate no. of days exceeding threshold if required
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
                        
                        %                         % Create empty array for output
                        %                         data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.AbsThresh));
                        
                        % Calculate number of days exceeding threshold
                        disp(['Calculating number of days when ',Variable,' exceeds ',num2str(inputs.AbsThresh(EA_3))])
                        data_calc(:,:,s,n_output) = nansum(data>inputs.AbsThresh(EA_3),3)/tslength;
                        titles2(n_output) = {[Variable, '>',num2str(inputs.AbsThresh(EA_3))]};
                        titles1(n_output) = {Dataset};
                        EA_3 = EA_3 + 1;
                    end
                            
                            
                    
%                     %% Calculate extreme mean if required
%                     if isfield(inputs,'Pctile')
%                         
% %                         % Create empty array for output
% %                         data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.Pctile));
%                         
%                         % Calculate given percentile
%                         for p = 1:length(inputs.Pctile)
%                             disp(['Calculating ',num2str(inputs.Pctile(p)),'th percentile over time period'])
%                             data_calc(:,:,s,n_output) = prctile(data,inputs.Pctile(p),3);
%                             titles2(p) = {[Variable, ' ',num2str(inputs.Pctile(p)),'th percentile']};
%                             titles1(p) = {Dataset};
%                             n_output = n_output + 1;
%                         end
%                     end
                            
%                     %% Calculate extreme mean if required
%                     if isfield(inputs,'ExtremeMeanPctile')
%                         
% %                         % Create empty array for output
% %                         data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.ExtremeMeanPctile));
%                         
%                         for p = 1:length(inputs.ExtremeMeanPctile)
%                             disp(['Calculating extreme mean above ',num2str(inputs.ExtremeMeanPctile(p)),'th percentile over time period'])
%                             
%                             % Calculate percentile threshold
%                             Txx_temp = data >= prctile(data,inputs.ExtremeMeanPctile(p),3);
%                             Txx = nan(size(Txx_temp));
%                             Txx(Txx_temp == 1) = 1;
%                             data_Txx = data .* Txx;
%                             % Find mean of days exceeding the threshold
%                             data_calc(:,:,s,n_output) = squeeze(nanmean(data_Txx,3));
%                             titles2(p) = {[Variable, ' extreme mean >',num2str(inputs.ExtremeMeanPctile(p)),'th percentile']};
%                             titles1(p) = {Dataset};
%                             n_output = n_output + 1;
%                         end
%                     end
                    
%                     %% Calculate no. of days exceeding threshold if required
%                     if isfield(inputs,'AbsThresh')
%                         
%                         % Find length of loaded time series (assuming whole summers
%                         % have been taken in inputs.TemporalRange)
%                         startyr = num2str(inputs.TemporalRange(1));
%                         if str2double(startyr(5:6))<=6
%                             startyr = str2double(startyr(1:4))-1;
%                         else
%                             startyr = str2double(startyr(1:4));
%                         end
%                         
%                         endyr = num2str(inputs.TemporalRange(2));
%                         if str2double(endyr(5:6))>=7
%                             endyr = str2double(endyr(1:4));
%                         else
%                             endyr = str2double(endyr(1:4))-1;
%                         end
%                         
%                         tslength = endyr-startyr;
%                         
% %                         % Create empty array for output
% %                         data_calc = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs.AbsThresh));
%                         
%                         % Calculate number of days exceeding threshold
%                         for p = 1:length(inputs.AbsThresh)
%                             disp(['Calculating number of days when ',Variable,' exceeds ',num2str(inputs.AbsThresh(p))])
%                             data_calc(:,:,s,n_output) = nansum(data>inputs.AbsThresh(p),3)/tslength;
%                             titles2(p) = {[Variable, '>',num2str(inputs.AbsThresh(p))]};
%                             titles1(p) = {Dataset};
%                             n_output = n_output + 1;
%                         end
%                     end
                    
                    
                    %% Plotting
                    % If calculating MMM or MMP, only plot all ensemble members if requested
                    if inputs.PlotAll ~= 0
%                         for p = 1:length(data_calc(1,1,:))
                            figure
                            UK_subplot(data_calc(:,:,s,n_output) .* LSM(grid_idx,grid_idy),[titles1(n_output),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
%                         end
                    end
                    
                    
                    % Calculate MMM or MMP
                    % Only do this once all data has been loaded
                    if s == length(inputs.Dataset)
                        
                        if inputs.MMM == 1
                            % Only want to mean model simulations, not obs. etc.
                            if sum(MM_id == 1)>1
                                data_plot = nanmean(data_calc(:,:,MM_id == 1,n_output),3);
                            elseif sum(MM_id == 2)>1
                                data_plot = nanmean(data_calc(:,:,MM_id == 2,n_output),3);                                
                            elseif sum(MM_id == 3)>1
                                data_plot = nanmean(data_calc(:,:,MM_id == 3,n_output),3);
                            end
                            
                            title1 = {'MMM'};
                            
                            % Plot
                            figure
                            UK_subplot(data_plot .* LSM(grid_idx,grid_idy),[title1(1),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                            
                        end
                            
                        if ~isempty(inputs.MMP)
                            
                            % Calculate MMP
                            for P = 1:length(inputs.MMP)
                                % Only want to mean model simulations, not obs. etc.
                                if sum(MM_id == 1)>1
                                    data_plot = prctile(data_calc(:,:,MM_id == 1,n_output),inputs.MMP(P),3);
                                elseif sum(MM_id == 2)>1
                                    data_plot = prctile(data_calc(:,:,MM_id == 2,n_output),inputs.MMP(P),3);
                                elseif sum(MM_id == 3)>1
                                    data_plot = prctile(data_calc(:,:,MM_id == 3,n_output),inputs.MMP(P),3);
                                end
                                
                                title1 = {['MM ',num2str(inputs.MMP(P)),'th percentile']};
                                
                                % Plot
                                figure
                                UK_subplot(data_plot .* LSM(grid_idx,grid_idy),[title1(1),titles2(n_output)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
                            
                            end
                        end

                    end
                        
                        
                        
                            
                            
%                             data_mod = data_calc;
%                         else
%                             data_mod = nan(size(data_calc));
%                         end
%                         
%                         
%                         % If this is the first simulation loaded, make new array
%                         if s == 1
%                             data_all = data_mod;
%                             
%                             % Otherwise, add the current simulation to the array
%                         elseif s > 1 && s < length(inputs.Dataset)
%                             data_all = cat(ndims(data_calc)+1,data_all,data_mod);
%                             
%                             % And if this is the final simulation, calculate the MMM or MMP
%                         elseif s == length(inputs.Dataset)
%                             data_all = cat(ndims(data_calc)+1,data_all,data_mod);
%                             
%                             % Calculate MMM if required
%                             if inputs.MMM == 1
%                                 data_plot = nanmean(data_all,ndims(data_calc)+1);
%                                 title1 = {'MMM'};
%                             end
%                                 
%                             % Calculate MMP if required    
%                             if ~isempty(inputs.MMP)
%                                 
%                                 % Create empty array for each MMP, if MMM not already calculated
%                                 if ~exist('data_plot','var')
%                                     data_plot = nan([size(data_calc),length(inputs.MMP)]);
%                                 % If MMM is already calculcated, add this to the end of the final plotting array    
%                                 else
%                                     data_plot_temp = nan([size(data_calc),length(inputs.MMP)+1]);
%                                     data_plot_temp(:,:,:,length(inputs.MMP)+1) = data_plot;
%                                     data_plot = data_plot_temp;
%                                     title1 = cell(1,length(inputs.MMP)+1);
%                                     title1(length(inputs.MMP)+1) = {'MMM'};
%                                 end
%                                 
%                                 % Calculate MMP
%                                 for P = 1:length(inputs.MMP)
%                                     data_plot(:,:,:,P) = prctile(data_all,inputs.MMP(P),ndims(data_calc)+1);
%                                     title1(P) = {['MM ',num2str(inputs.MMP(P)),'th percentile']};
%                                 end
%                             end
%                             
%                             % Go through each precentile to plot MMM/MMP
%                             for p = 1:length(data_calc(1,1,:))
%                                 for P = 1:length(data_plot(1,1,1,:))
%                                     figure
%                                     UK_subplot(data_plot(:,:,p,P) .* LSM(grid_idx,grid_idy),[title1(P),titles2(p)],[],lats(grid_idx,grid_idy),lons(grid_idx,grid_idy),inputs)
%                                 end
%                             end
%                         end
%                     end
                end
            end
        end
    end
end


disp('-----')


