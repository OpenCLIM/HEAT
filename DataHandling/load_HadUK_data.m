function [data80s,data90s,data00s,data10s] = load_HadUK_data(var,res,summer,specdecade)
% [data80s,data90s,data00s,data10s] = load_HadUK_data(var,res,summer,specdecade)
%
% Load the HadUK-Grid observational data for UK. The 1km dataset is too large
% to save the entire 39 year climatology (even just summers only) to one
% matrix, so it is loaded in 10 (or 8) year chunks. 
% 
%
% Output:
%   data80s = the loaded data for 1981-1990
%   data90s = the loaded data for 1991-2000
%   data00s = the loaded data for 2001-2010
%   data10s = the loaded data for 2011-2018 (end of available data)
%
% Input:
%   var = variable to be loaded ('tas', 'tasmax', 'tasmin' 'pv', 'latitude',' longitude' or 'time')
%   res = the resolution to be loaded ('1km', '12km' or '60km')
%   summer = 'JJA' or 'MO' [default] (Nb. annual is not an option for all
%       varaibles as I did not download all months)
%   specdecade = if only want to load a specific decade, define it here as
%       '1980s',' 1990s', '2000s' or 2010s', then data80s output var will
%       be only this decade (optional)
%

%% Set defaults
if ~exist('summer','var')
    summer = 'MO';
end


%% Find which machine is being used
curdir = pwd;

% Run locally
if strcmp(curdir(1:14),'/Users/ak0920/') % Currently just works locally
    disp('load_HadUK_data.m: Running locally')
    % Set data dir
    data_dir = '/Volumes/DataDrive/HadUK-Grid/v1.0.2.1/';
    
else % Run on Anthropocene
    if srtcmp(curdir(1:14),'/home/bridge/a')
        data_dir = '/export/anthropocene/array-01/ak0920/HadUKGrid/';
    end
end


%% Generate filename that is to be loaded
% Some variables (lat, lon and time) need loaded from a different variable file
if strcmp(var,'latitude') || strcmp(var,'longitude') || strcmp(var,'time')
    var2 = 'tasmax';
    catdir = 1;
else
    var2 = var;
    catdir = 3;
end

% Set resolution and variable specific file directory
filedir = [var2,'/',res,'/'];

% Define the years for each decade
years = [1981:1990;1991:2000;2001:2010;2011:2019,nan];

% If only interested in a specific decade
if exist('specdecade','var')
    if strcmp(specdecade,'1980s')
        years = years(1,:);
    else
        if strcmp(specdecade,'1990s')
            years = years(2,:);
        else
            if strcmp(specdecade,'2000s')
                years = years(3,:);
            else
                if strcmp(specdecade,'2010s')
                    years = years(4,:);
                end
            end
        end
    end
end

%% Load data
% if strcmp(res,'1km')

% Go through each decade
% for decade = 1:length(years(:,1))
for decade = 1:nargout
    
    % Go through years of each decade
    for year = years(decade,:)
        
        % Use this later to figure out if year is first of decade
        year1 = num2str(year);
        
        if ~isnan(year)
            % Sanity check:
            disp(['Loading year: ',num2str(year)])
            
            %% Load tasmin or tasmax
            % These have different file structure as they are saved monthly, not for whole year:
            if strcmp(var2,'tasmax') || strcmp(var2,'tasmin')
                
                % Daily temperature range file
                file_june = [filedir,var2,'_hadukgrid_uk_',res,'_day_',num2str(year),'0601-',num2str(year),'0630.nc'];
                file_july = [filedir,var2,'_hadukgrid_uk_',res,'_day_',num2str(year),'0701-',num2str(year),'0731.nc'];
                file_aug = [filedir,var2,'_hadukgrid_uk_',res,'_day_',num2str(year),'0801-',num2str(year),'0831.nc'];
                file_sept = [filedir,var2,'_hadukgrid_uk_',res,'_day_',num2str(year),'0901-',num2str(year),'0930.nc'];
                
                % Load data
                data_june = ncread([data_dir,file_june],var);
                data_july = ncread([data_dir,file_july],var);
                data_aug = ncread([data_dir,file_aug],var);
                data_sept = ncread([data_dir,file_sept],var);
                
                
                % Time data is only in 1D - select only the 15 days of
                % September that are needed:
                if strcmp(var,'time')
                    data_sept = data_sept(1:15);
                else
                    data_sept = data_sept(:,:,1:15);
                end
                
                
                % Concatenate into one file for correct length of summer
                if strcmp(summer,'MO')
                    if strcmp(year1(4),'1') % if first year in decade start new 'data'
                        data = cat(catdir,data_june,data_july,data_aug,data_sept);
                    else % otherwise add year onto existing 'data'
                        data = cat(catdir,data,data_june,data_july,data_aug,data_sept);
                    end
                else
                    if strcmp(summer,'JJA')
                        if strcmp(year1(4),'1')
                            data = cat(catdir,data_june,data_july,data_aug);
                        else
                            data = cat(catdir,data,data_june,data_july,data_aug);
                        end
                    end
                end
                
                
            else
                %% Load any other variable (saved annually)
                file = [filedir,var2,'_hadukgrid_uk_',res,'_mon_',num2str(year),'01-',num2str(year),'12.nc'];
                data_ann = ncread([data_dir,file],var2);
                
                % Save/concatenate monthly as daily data
                if strcmp(summer,'MO')
                    if strcmp(year1(4),'1') % if first year in decade start new 'data'
                        data = cat(catdir,repmat(data_ann(:,:,6),1,1,30),repmat(data_ann(:,:,7),1,1,31),repmat(data_ann(:,:,8),1,1,31),repmat(data_ann(:,:,9),1,1,15));
                    else % otherwise add year onto existing 'data'
                        data = cat(catdir,data,repmat(data_ann(:,:,6),1,1,30),repmat(data_ann(:,:,7),1,1,31),repmat(data_ann(:,:,8),1,1,31),repmat(data_ann(:,:,9),1,1,15));
                    end
                else
                    if strcmp(summer,'JJA')
                        if strcmp(year1(4),'1') % if first year in decade start new 'data'
                            data = cat(catdir,repmat(data_ann(:,:,6),1,1,30),repmat(data_ann(:,:,7),1,1,31),repmat(data_ann(:,:,8),1,1,31));
                        else % otherwise add year onto existing 'data'
                            data = cat(catdir,data,repmat(data_ann(:,:,6),1,1,30),repmat(data_ann(:,:,7),1,1,31),repmat(data_ann(:,:,8),1,1,31));
                        end
                    end
                end
                
            end
            
        end
    end
    
    
    %% Tidy up and correct units
    if ~strcmp(var,'time')
        data(data>1000) = nan;
    end
    
    %% Save data to the correct decade
    if decade == 1
        data80s = data;
    else
        if decade == 2
            data90s = data;
        else
            if decade == 3
                data00s = data;
            else
                data10s = data;
            end
        end
    end
    
end


