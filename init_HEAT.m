% init_HEAT.m
%
% Intialise the HEAT tool box. Sets up directory paths to core climate
% input data and adds any necessary paths.

%% Find which machine is being used
root_dir = pwd;

% Running locally on AKA's machine
if strcmp(root_dir(1:14),'/Users/ak0920/')
%     disp('Initialising HEAT on local machine')
    
    % Set data directories
    UKCP18dir = '/Volumes/DataDrive/UKCP18/';
    ERA5dir = '/Volumes/DataDrive/ERA5/';
    HadUKdir = '/Volumes/DataDrive/HadUK-Grid/';
    Deriveddir = '/Volumes/DataDrive/HEAToutput/';
    addpath('PhysicalCalculations/')
    addpath('DataHandling/')
    addpath('Processing/')
    addpath('Outputting/')
    addpath('PreProcessedData/')
    
    % Running on Anthropocene
else
    if strcmp(root_dir(1:14),'/home/bridge/a')
%         disp('Initialising HEAT on Anthropocene')
        
        % Set data directory
        data_dir = '/export/anthropocene/array-01/ak0920/ukcp18_data/';
        
        % Running on DAFNI
    else
        if strcmp(root_dir(1:14),'/DAFNI/dir/TBC/')
%             disp('Initialising HEAT on DAFNI')
            
            % Set data directory
            data_dir = '/DAFNI/data/dir/TBC/';
            
            % Running on another machine
        else
            disp('Initialising HEAT on unknown machine')
            prompt = 'Please enter UKCP18 root data directory: ';
            data_dir = input(prompt,'s');
        end
    end
end

