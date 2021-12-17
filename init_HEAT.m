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
        addpath('Dependencies/')
        
        
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
        addpath('Dependencies/')
        
        
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
        addpath('Dependencies/')
        
    elseif strcmp(root_dir(1:14),'/home/alanka/H')
        UKCP18dir = '/data/ukcp18_data/';
        ERA5dir = '/data/ERA5/';
        HadUKdir = '/data/HadUK-Grid/';
        Deriveddir = '/data/HEAToutput/DerivedData/';
        Outputdir = '/data/HEAToutput/';
        addpath('PhysicalCalculations/')
        addpath('DataHandling/')
        addpath('Processing/')
        addpath('Outputting/')
        addpath('PreProcessedData/')
        addpath('Dependencies/')
        
    end
    % Otherwise, assume running on DAFNI
else
    % Set data directory
    cd('/data/')
    UKCP18dir = '/inputs/UKCP18dir/';
    Outputdir = '/outputs/';
    
    % ATK-A I hope that these directories can be mounted to using:
    %     docker run -v /Path/to/UKCP18/on/DAFNI/:/data/UKCP18dir/ -v /Path/to/output/storage/on/DAFNI/:/data/HEAToutput/ heat
    
    
    % Running on another machine
    
    
%     addpath('PhysicalCalculations/')
%     addpath('DataHandling/')
%     addpath('Processing/')
%     addpath('Outputting/')
%     addpath('PreProcessedData/')
%     addpath('Dependencies/')
    
    % Note: This needs updated so that the entered directories become
    % global and do not need re-entered everytime init_HEAT.m is called
    % (which is quite often).
    
end
