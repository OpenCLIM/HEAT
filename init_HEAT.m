% init_HEAT.m
%
% Intialise the HEAT tool box. Sets up directory paths to core climate
% input data and adds any necessary paths.

% Find which machine is being used
root_dir = pwd;

% Running locally on AKA's machine
if length(root_dir) >= 14
    if strcmp(root_dir(1:14),'/Users/ak0920/')
        % Set data directories
%         UKCP18dir = '/Volumes/DataDrive/UKCP18/';
%         UKCP18dir = '/Users/ak0920/Data/';
%         Climatedirin = '/Volumes/BCDATA/';
%         ERA5dir = '/Volumes/DataDrive/ERA5/';
%         HadUKdir = '/Volumes/DataDrive/HadUK-Grid/v1.0.2.1/';
        Deriveddir = '/Volumes/DataDrive/HEAToutput/DerivedData/';
        Inputdir = '/data/inputs/';
        Climatedirout = '/Volumes/DataDrive/HEAToutput/';
%         Climatedirout = '/Users/ak0920/Data/HEAToutput/';
        
        addpath('PhysicalCalculations/')
        addpath('DataHandling/')
        addpath('Processing/')
        addpath('Outputting/')
        addpath('PreProcessedData/')
        addpath('Dependencies/')
        
    elseif strcmp(root_dir(1:14),'/home/alanka/H')
        Climatedirin = '/data/ukcp18_data/';
        ERA5dir = '/data/ERA5/';
        HadUKdir = '/data/HadUK-Grid/';
        Deriveddir = '/data/HEAToutput/DerivedData/';
        Climatedirout = '/data/HEAToutput/';
        addpath('PhysicalCalculations/')
        addpath('DataHandling/')
        addpath('Processing/')
        addpath('Outputting/')
        addpath('PreProcessedData/')
        addpath('Dependencies/')
        
    end
    
    
    % Otherwise, assume running on DAFNI (in a Docker container)
else
    % Set input and output data directories
    Inputdir = '/data/inputs/';
    Climatedirin = '/data/inputs/ClimateData/';
    Urbandirin = '/data/inputs/Urban/';
    SubsetClimatedirin = '/data/inputs/ClimateAlt/';
    BaseClimatedirin = '/data/inputs/ClimateBase/';
    
    Climatedirout = '/data/outputs/Climate/';


    % Note: in this case, 'addpath' should have been done before building
    % the Matlab executable and Docker file, by running add_paths.m.
    
end
