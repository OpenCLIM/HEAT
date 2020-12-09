function [data,yyyymmdd] = load_UKCP_data(model,var,CPM_period,reg)
% [data] = load_UKCP_data(model,var,CPM_period,reg)
%
% This function loads UKCP18 data, specified by simulation, variable,
% temporal frequency and spatial coverage.
%
% Outputs:
%   data = the desired dataset
%   yyyymmdd = date information to allow correction of some missing date
%       values in UKCP18 CMIP5 member 16
%
% Inputs:
%   model = string, 'gcm...', 'rcm...' or 'cpm...' where '...' is the scenario (e.g.
%       85 for RCP 8.5 follwed by ensemble member (01-28 for GCMs, 01-15
%       for RCMs), e.g. 'gcm8501'
%   var = string, variable to be loaded in original naming convention:
%       T = 'tasmax', q = 'huss', p = 'psl', time step date = 'yyyymmdd'
%   reg = spatial domaing: 'uk' [default] or 'global'
%   CPM_period = the start of the time period of CPM data to load
%       (1981 [1981-2000], 2021 [2012-2040], 2061 [2061-20-80])



%% Set some defaults
%  CPM time period
if ~exist('CPM_period','var')
    CPM_period = 1981';
end

% Domain to load
if ~exist('reg','var')
    reg = 'uk';
end

% Temporal frequency to load (currently always daily)
if ~exist('tempfreq','var')
    tempfreq = 'day';
end

% Convert resolution description for file directory description
if strcmp(model(1:3),'gcm')
    res = '60km';
else
    if strcmp(model(1:3),'rcm')
        res = '12km';
    else
        res = '2km';
    end
end


%% Corrections for loading time step/date info
% If reading 'yyyymmdd' or 'time', it should be taken from a netCDF file
% within the folder for a different variable (i.e. yyyymmdd is saved within
% a tas file, and does not have a directory path and file dedicated to it).
% Also, as time step is 1D, the data should be compiled from multiple files
% in 2nd, not 3rd, dimension
if strcmp(var,'yyyymmdd')
    var1 = 'tas';
    catdim = 2;
else
    if strcmp(var,'time') && strcmp(reg,'global') % time variables have different diemsion structure for daily and monthly files
        var1 = 'tas';
        catdim = 2;
    else
        if strcmp(var,'time') && strcmp(reg,'uk') % time variables have different diemsion structure for daily and monthly files
            var1 = 'tas';
            catdim = 1;
        else % all other variables are fine
            var1 = var;
            catdim = 3;
        end
        
    end
end

% If loading time step with monthly data, the variable name should be
% 'yyyymm', not 'yyyymmdd': correct this if necessary
if strcmp(var,'yyyymmdd') && strcmp(tempfreq,'mon')
    var = 'yyyymm';
end

%% Find which machine is being used
curdir = pwd;

%% Run locally
if strcmp(curdir(1:14),'/Users/ak0920/')
%     disp(['load_UKCP_data.m: ',model])
    
    % Find data directory
    if strcmp(reg,'global')
        data_dir = '/Volumes/DataDrive/UKCP18/tas_monthly/';
    else
        data_dir = ['/Volumes/DataDrive/UKCP18/',res,'/',var1,'/run',model(6:7),'/'];
    end
    
    % Days to load from cpm data as loading whole year makes file too big
    if strcmp(var,var1)
        starts = [1 1 121 1]; % April to September
        ends = [Inf Inf 180 Inf];
    else
        if strcmp(var,'yyyymmdd')
            starts = [1 121];
            ends = [Inf 180];
        else
            starts = 121;
            ends = 180;
        end
    end
    
    % Load CPM data
    if strcmp(model(1:3),'cpm')
        time_starts = CPM_period-1:CPM_period+18;
        
        % Load each time slice
        for i = 1:length(time_starts)
            
            disp(['Year: ',num2str(time_starts(i)+1)])
            
            filename = [data_dir,var1,'_rcp',model(4:5),'_land-cpm_uk_2.2km_',model(6:7),...
                '_',tempfreq,'_',num2str(time_starts(i)),'1201-',num2str(time_starts(i)+1),'1130.nc'];
            
            % Load data and concatenate into one array
            if i == 1
                data = ncread(filename,var,starts,ends);
                yyyymmdd = ncread(filename,'yyyymmdd',[1 121],[Inf 180]);
            else
                data = cat(catdim,data,ncread(filename,var,starts,ends));
                yyyymmdd = cat(2,yyyymmdd,ncread(filename,'yyyymmdd',[1 121],[Inf 180]));
            end
        end
    else
        
        
        % Load RCM data
        if strcmp(model(1:3),'rcm')
            
            % Other variables are in decadal time slices
            time_starts = 1980:10:2070;
            
            % Load each time slice
            for i = 1:length(time_starts)
                
                filename = [data_dir,var1,'_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),...
                    '_',tempfreq,'_',num2str(time_starts(i)),'1201-',num2str(time_starts(i)+10),'1130.nc'];
                
                % Load data and concatenate into one array
                if i == 1
                    data = ncread(filename,var);
                    yyyymmdd = ncread(filename,'yyyymmdd');
                else
                    data = cat(catdim,data,ncread(filename,var));
                    yyyymmdd = cat(2,yyyymmdd,ncread(filename,'yyyymmdd'));
                end
            end
            
            % Load GCM data
        else
            % Time period of netCDF depends on whether monthly or daily is used
            if strcmp(tempfreq,'day')
                time_starts = 1979:10:2079;
                
                for i = 1:length(time_starts)
                    filename = [data_dir,var1,'_rcp',model(4:5),'_land-gcm_',reg,'_60km_',model(6:7),...
                        '_',tempfreq,'_',num2str(time_starts(i)),'1201-',num2str(time_starts(i)+10),'1130.nc'];
                    
                    % Load data and concatenate into one array
                    if i == 1
                        data = ncread(filename,var);
                        yyyymmdd = ncread(filename,'yyyymmdd');
                    else
                        data = cat(catdim,data,ncread(filename,var));
                        yyyymmdd = cat(2,yyyymmdd,ncread(filename,'yyyymmdd'));
                    end
                end
            else % Monthly data is all in one file
                filename = [data_dir,var1,'_rcp',model(4:5),'_land-gcm_',reg,'_60km_',model(6:7),...
                    '_',tempfreq,'_189912-209911.nc'];
                % Load the data (concatenation not needed)
                data = ncread(filename,var);
            end
        end
        
    end
    
    
%% Run on Anthropocene
else
    if srtcmp(curdir(1:14),'/home/bridge/a')
%         disp(['load_UKCP_data.m: ',model])
        
        % Find data directory:
        % Some variables are stored in my directory, others in Eunice's.
        % The psl variable for the GCM has a different directory path structure in
        % Eunice's directory and her tas and tasmax netCDF files have been merged
        % from decadal into one long time series. This must all be accounted for...
        
        if strcmp(model(1:3),'rcm')
            if strcmp(var,'huss')
                data_dir = ['/export/anthropocene/array-01/ak0920/ukcp18_data/12km/',...
                    var1,'/',tempfreq,'/run',model(6:7),'/'];
                files_merged = 0;
            else
                data_dir = ['/export/anthropocene/array-01/yl17544/ukcp18_data/12km/',...
                    var1,'/',tempfreq,'/run',model(6:7),'/'];
                if strcmp(var,'psl')
                    files_merged = 0;
                else
                    files_merged = 1; % Eunice has merged the tas and tasmax RCM netCDF files, I have not
                end
            end
        else
            if strcmp(var,'psl')
                data_dir = ['/export/anthropocene/array-01/yl17544/ukcp18_data/60km/',tempfreq,'/run',...
                    model(6:7),'/'];
                files_merged = 0;
            else
                data_dir = ['/export/anthropocene/array-01/ak0920/ukcp18_data/60km/',...
                    var1,'/',tempfreq,'/run',model(6:7),'/'];
                files_merged = 0;
            end
        end
        
        
        % Load tas or tasmax RCM data from Eunice's directory
        % (i.e. those that are already merged)
        if files_merged == 1
            
            filename = [data_dir,var1,'_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),...
                '_',tempfreq,'_19801201-20801130.nc'];
            data = ncread(filename,var);
            
            
            % Otherwise load the other variables
        else
            % Load RCM data
            if strcmp(model(1:3),'rcm')
                
                % Other variables are in decadal time slices
                time_starts = 1980:10:2070;
                
                % Load each time slice
                for i = 1:length(time_starts)
                    
                    filename = [data_dir,var1,'_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),...
                        '_',tempfreq,'_',num2str(time_starts(i)),'1201-',num2str(time_starts(i)+10),'1130.nc'];
                    
                    % Load data and concatenate into one array
                    if i == 1
                        data = ncread(filename,var);
                    else
                        data = cat(catdim,data,ncread(filename,var));
                    end
                end
                
                % Load GCM data
            else
                time_starts = 1979:10:2079;
                
                for i = 1:length(time_starts)
                    filename = [data_dir,var1,'_rcp',model(4:5),'_land-gcm_',reg,'_60km_',model(6:7),...
                        '_',tempfreq,'_',num2str(time_starts(i)),'1201-',num2str(time_starts(i)+10),'1130.nc'];
                    
                    % Load data and concatenate into one array
                    if i == 1
                        data = ncread(filename,var);
                    else
                        data = cat(catdim,data,ncread(filename,var));
                    end
                end
            end
        end
    end
end

%% Convert pressure units in RCM simulations
if strcmp(model(1),'r') && strcmp(var,'psl')
    data = data/100;
end

end


