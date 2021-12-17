%% Description
% This .m file provides the input file for a default test of the HEAT code
% when importing to DAFNI. It produces some simple output including a map
% and ARCADIA-ready .csvs using the GCMs from UKCP18.
% 
% The output should consist of:
%  1) Maps of number of days exceeding Tmean>25, one for each of
%   the two models and one multi-model mean. This will be saved in
%   HEAToutput/HEAT_DAFNI_test/.
%  2) 391 .csv files, one for each cell of the GCM grid, each with 15 rows
%   (one for each model) and 10800 columns (daily data for 30 years with 360
%   day calendar year). This will be saved in HEAToutput/HEAT_DAFNI_test/.
% 

%% Set some basic experiment details
inputs.ExptName = 'HEAT_DAFNI_test'; 
inputs.OverwriteExpt = 1; % OR inputs.OverwriteExpt = 1;
inputs.Domain = 'UK'; % OR inputs.Domain = 'global';
inputs.SaveDerivedOutput = 0; % OR inputs.SaveDerivedOutput = 0;
inputs.OverwriteDerivedOutput = 0; % OR inputs.OverwriteDerivedOutput = 1;
% inputs.BiasCorr = 0; % OR inputs.BiasCorr = 0; 


%% Select dataset(s) to use
inputs.DataType = {'UKCP18'};
inputs.Variable = {'Tmean'};
% inputs.Dataset = {'GCM-01','GCM-02','GCM-03','GCM-04','GCM-05','GCM-06','GCM-07','GCM-08','GCM-09','GCM-10','GCM-11','GCM-12','GCM-13','GCM-14','GCM-15'};
inputs.Dataset = {'GCM-01'};
% inputs.MMM = 1;


%% Subsetting of data for analysis
inputs.TemporalRange = [20500101, 20791230]; % yyyy, yyyymm or yyyymmdd start and end dates
% inputs.AnnSummer = 'Summer'; % OR inputs.AnnSummer = 'Summer'; 'JJA'; 'MJJAS','ann'; % Summer = 1st June-15th Sept. or leave blank for annual mean
% inputs.SpatialRange = [51,55.5;-11, -5]; % [latS,latN;lonW,lonE] for boxed region or [lat,lon] for single grid point
% inputs.Region = {}; % 'Scotland','North East','North West','Yorkshire and 
%     % the Humber','East Midlands','West Midlands','East of England','Greater
%     % London','South East','South West','Wales','Northern Ireland','Isle of
%     % Man','RoI','Ireland','GB'


%% Output types
% inputs.AbsThresh = 25;
% inputs.OutputType = {'map'};
% inputs.AveTime = 10; % Time series averaging length in years: default = 10
inputs.SaveFigures = 1;
inputs.PlotAll = 1;
% inputs.OutputRegion = 'all';
inputs.WorkflowOutput = 'ARCADIA';

