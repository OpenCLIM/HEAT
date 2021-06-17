%% Description
% This .m file provides the input file for a default test of the HEAT code
% when importing to a new machine. It produces a simple output from each
% step of the code (a derived netCDF, a map and ARCADIA-ready .csvs) using
% two GCMs from UKCP18.
% 
% The output should consist of:
%  1) New derived data netCDF files for sWBGTmax for GCM-01 and GCM-04 out
%   of the UKCP18 ensemble, for as many years as there is data available.
%   This will be saved in HEAToutput/DerivedData/.
%  2) Three maps of number of days exceeding sWBGTmax>25, one for each of
%   the two models and one multi-model mean. This will be saved in
%   HEAToutput/HEAT_test/.
%  3) 391 .csv files, one for each cell of the GCM grid, each with 2 rows
%   (one for each model) and 3600 columns (daily data for 10 years with 360
%   day calendar year). This will be saved in HEAToutput/HEAT_test/.
% 

%% Set some basic experiment details
inputs.ExptName = 'HEAT_test'; 
inputs.OverwriteExpt = 1; % OR inputs.OverwriteExpt = 1;
inputs.Domain = 'UK'; % OR inputs.Domain = 'global';
inputs.SaveDerivedOutput = 1; % OR inputs.SaveDerivedOutput = 0;
inputs.OverwriteDerivedOutput = 1; % OR inputs.OverwriteDerivedOutput = 1;
% inputs.BiasCorr = 0; % OR inputs.BiasCorr = 0; 


%% Select dataset(s) to use
inputs.DataType = {'UKCP18'};
inputs.Variable = {'sWBGTmax'};
inputs.Dataset = {'GCM-01','GCM-04'};
inputs.MMM = 1;


%% Subsetting of data for analysis
inputs.TemporalRange = [20060101, 20151230]; % yyyy, yyyymm or yyyymmdd start and end dates
inputs.AnnSummer = 'Summer'; % OR inputs.AnnSummer = 'Summer'; 'JJA'; 'MJJAS','ann'; % Summer = 1st June-15th Sept. or leave blank for annual mean
% inputs.SpatialRange = [51,55.5;-11, -5]; % [latS,latN;lonW,lonE] for boxed region or [lat,lon] for single grid point
% inputs.Region = {}; % 'Scotland','North East','North West','Yorkshire and 
%     % the Humber','East Midlands','West Midlands','East of England','Greater
%     % London','South East','South West','Wales','Northern Ireland','Isle of
%     % Man','RoI','Ireland','GB'


%% Output types
inputs.AbsThresh = 25;


inputs.OutputType = {'map'};
% inputs.AveTime = 10; % Time series averaging length in years: default = 10
inputs.SaveFigures = 1;
inputs.PlotAll = 1;
% inputs.OutputRegion = 'all';
inputs.WorkflowOutput = 'ARCADIA';

