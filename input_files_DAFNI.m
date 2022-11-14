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
inputs.ExptName = 'HEAT_DAFNI_default';
inputs.OverwriteExpt = 1; % OR inputs.OverwriteExpt = 1;
inputs.Domain = 'UK'; % OR inputs.Domain = 'global';
inputs.SaveDerivedOutput = 0; % OR inputs.SaveDerivedOutput = 0;
inputs.OverwriteDerivedOutput = 0; % OR inputs.OverwriteDerivedOutput = 1;
inputs.BiasCorr = 0; % OR inputs.BiasCorr = 0;


%% Select dataset(s) to use
inputs.DataType = {'UKCP18'};
inputs.Variable = {'Tmean'};
inputs.Variable = {'tas'};
inputs.Dataset = {'RCM-01','RCM-04','RCM-05','RCM-06','RCM-07','RCM-08','RCM-09','RCM-10','RCM-11','RCM-12','RCM-13','RCM-15'};
% inputs.Dataset = {'GCM-01','GCM-02','GCM-03','GCM-04','GCM-05','GCM-06','GCM-07','GCM-08','GCM-09','GCM-10','GCM-11','GCM-12','GCM-13','GCM-14','GCM-15'};
% inputs.MMM = 1;


%% Subsetting of data for analysis
inputs.Scenario = 'past';
% inputs.TemporalRange = [20100101, 20181230]; % yyyy, yyyymm or yyyymmdd start and end dates
inputs.PeriodStart = 1990; % yyyy, yyyymm or yyyymmdd start and end dates
% inputs.AnnSummer = 'Summer'; % OR inputs.AnnSummer = 'Summer'; 'JJA'; 'MJJAS','ann'; % Summer = 1st June-15th Sept. or leave blank for annual mean
% inputs.SpatialRange = [51,55.5;-11, -5]; % [latS,latN;lonW,lonE] for boxed region or [lat,lon] for single grid point
% inputs.Region = {}; % 'Scotland','North East','North West','Yorkshire and
%     % the Humber','East Midlands','West Midlands','East of England','Greater
%     % London','South East','South West','Wales','Northern Ireland','Isle of
%     % Man','RoI','Ireland','GB'


%% Parameters
inputs.MMTpctile = 93;
% inputs.AbsThresh = 25;
% inputs.AveTime = 10; % Time series averaging length in years: default = 10
% inputs.ExtremeMeanPctile = [95 99];


%% Output types
inputs.SaveFigures = 1;
inputs.PlotAll = 1;
% inputs.OutputRegion = 'all';
% inputs.OutputType = 'NetCDF'; % {'map'};
% inputs.OutputType = 'Extreme mean'; % {'map'};
inputs.OutputType = 'DD66'; % {'map'};


%% Update fields if Environment variables are provided
disp('Updating inputs with environment variables')
% Then overwrite defaults with environment variables if running on DAFNI:
env_BC = getenv('BIASCORR');
env_expn = getenv('EXPNAME');
env_varn = getenv('VARNAME');
env_scen = getenv('SCENARIO');
env_tims = getenv('TIMEPERIOD_S');
env_timl = getenv('TIMEPERIOD_L');
env_outp = getenv('OUTPUT');
env_adap = getenv('ADAPT');
env_uhii = getenv('UHI_I');


if ~isempty(env_BC)
    if strcmp(env_BC,'y')
        disp('Environment variable found for bias correction option: updating inputs file')
        inputs.BiasCorr = 1;
    end
end

if ~isempty(env_expn)
    disp('Environment variable found for Experiment Name: updating inputs file')
    inputs.ExptName = {env_expn};
%     inputs.ExptName
end
if ~isempty(env_varn)
    disp('Environment variable found for Variable: updating inputs file')
    inputs.Variable = {env_varn};
%     inputs.Variable
end
if ~isempty(env_scen)
    disp('Environment variable found for Scenario: updating inputs file')
    inputs.Scenario = {string(env_scen)};
end

if ~isempty(env_tims)
    disp('Environment variable found for defining time period: updating inputs file')
    inputs.PeriodStart = env_tims;
end

if ~isempty(env_timl)
    inputs.PeriodLength = env_timl;
end

if ~isempty(env_outp)
    inputs.OutputType = env_outp;
end

if ~isempty(env_adap)
    inputs.MMTpctile = env_adap;
end

if ~isempty(env_uhii)
    inputs.UHI_I = env_uhii;
else
    inputs.UNI_I = 2;
end

disp(' ')

