% init_HEAT.m
%
% Intialise the HEAT tool box. Sets up directory paths to core climate
% input data and adds any necessary paths.

% Find which machine is being used
root_dir = pwd;

% Running locally on AKA's machine
if strcmp(root_dir(1:14),'/Users/ak0920/')
    % Set data directories
    UKCP18dir = '/Volumes/DataDrive/UKCP18/';
    ERA5dir = '/Volumes/DataDrive/ERA5/';
    HadUKdir = '/Volumes/DataDrive/HadUK-Grid/v1.0.2.1/';
    Deriveddir = '/Volumes/DataDrive/HEAToutput/DerivedData/';
    Outputdir = '/Volumes/DataDrive/HEAToutput/';
    addpath('PhysicalCalculations/')
    addpath('DataHandling/')
    addpath('Processing/')
    addpath('Outputting/')
    addpath('PreProcessedData/')
    
    
    % Running on Anthropocene
elseif strcmp(root_dir(1:14),'/home/bridge/a')
    % Set data directories
    UKCP18dir = '/export/anthropocene/array-01/ak0920/ukcp18_data/';
    ERA5dir = '/export/anthropocene/array-01/ak0920/ERA5/';
    HadUKdir = '/export/anthropocene/array-01/ak0920/HadUKGrid/';
    Deriveddir = '/export/anthropocene/array-01/ak0920/HEAToutput/DerivedData/';
    Outputdir = '/export/anthropocene/array-01/ak0920/HEAToutput/';
    addpath('PhysicalCalculations/')
    addpath('DataHandling/')
    addpath('Processing/')
    addpath('Outputting/')
    addpath('PreProcessedData/')
    
    
    % Running on BluePebble
elseif strcmp(root_dir(1:14),'/home/ak0920/h')
    % Set data directories
    UKCP18dir = '/bp1store/geog-tropical/data/ukcp18_data/';
    ERA5dir = '/bp1store/geog-tropical/data/ERA-5/';
    HadUKdir = '/bp1store/geog-tropical/data/HadUK-Grid/';
    Deriveddir = '/work/ak0920/HEAToutput/DerivedData/';
    Outputdir = '/work/ak0920/HEAToutput/';
    addpath('PhysicalCalculations/')
    addpath('DataHandling/')
    addpath('Processing/')
    addpath('Outputting/')
    addpath('PreProcessedData/')
    
    
    % Running on DAFNI
elseif strcmp(root_dir(1:14),'/DAFNI/dir/TBC/')
    % Set data directory
    data_dir = '/DAFNI/data/dir/TBC/';
    
    
    % Running on another machine
else
    
    disp('Initialising HEAT on unknown machine')
    
    prompt = 'Please enter UKCP18 root data directory: ';
    UKCP18dir = input(prompt,'s');
    prompt = 'Please enter ERA5 root data directory: ';
    ERA5dir = input(prompt,'s');
    prompt = 'Please enter HadUK-Grid root data directory: ';
    HadUKdir = input(prompt,'s');
    prompt = 'Please enter Derived data directory: ';
    Deriveddir = input(prompt,'s');
    prompt = 'Please enter Output data directory: ';
    Outputdir = input(prompt,'s');
    
    addpath('PhysicalCalculations/')
    addpath('DataHandling/')
    addpath('Processing/')
    addpath('Outputting/')
    addpath('PreProcessedData/')
    
end
