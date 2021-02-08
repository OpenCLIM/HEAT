%% Description
% This .m file is a template for all the required variables to lauch HEAT.
% A similar set of variables can be created using generate_launch_file.m
% using a simple GUI. 
% 
% This template allows variables to be selected without having to go
% through the GUI each time. However, certain combinations of variables
% cannot go together (e.g. there is no global humidity observational data)
% therefore these combinations cannot be loaded and processed. The GUI will
% only permit possible selections - using this template requires a careful
% choice of the correct variables.

% inputs1.ExptName = 'TemplateExpt'; % Give the experiment a name (no spaces)
% inputs1.Domain = 'UK'; % OR inputs1.Domain = 'global';
% inputs1.DataType = {'UKCP18','HadUKGrid','ERA5','CMIP6'};
% inputs1.Dataset = {'GCM-01','GCM-02','GCM-03','GCM-04','GCM-05','GCM-06','GCM-07','GCM-08','GCM-09','GCM-10','GCM-11','GCM-12','GCM-13','GCM-14','GCM-15',...
%             'CMIP5-16','CMIP5-19','CMIP5-21','CMIP5-23','CMIP5-24','CMIP5-25','CMIP5-27','CMIP5-28',...
%             'RCM-01','RCM-04','RCM-05','RCM-06','RCM-07','RCM-08','RCM-09','RCM-10','RCM-11','RCM-12','RCM-13','RCM-15',...
%             'CPM-01','CPM-04','CPM-05','CPM-06','CPM-07','CPM-08','CPM-09','CPM-10','CPM-11','CPM-12','CPM-13','CPM-15',...
%             '60km','12km','2km','1km','025deg'};
% inputs1.Variable = {'Tmean','Tmax','Tmin','VP','sWBGT','HD','AT'};
% inputs1.TempRes = 'daily'; % OR inputs1.TempRes = 'monthly';
% inputs1.BiasCorr = 1; % OR inputs1.BiasCorr = 0; 
% inputs1.SaveDerivedOutput = 1; % OR inputs1.SaveDerivedOutput = 0;
% inputs1.OverwriteDerivedOutput = 1; % OR inputs1.OverwriteDerivedOutput = 1;


%% Set some basic experiment details
inputs.ExptName = 'testing'; % Give the experiment a name (no spaces)
inputs.OverwriteExpt = 1; % OR inputs.OverwriteExpt = 1;
inputs.Domain = 'UK'; % OR inputs.Domain = 'global';
inputs.SaveDerivedOutput = 1; % OR inputs.SaveDerivedOutput = 0;
inputs.OverwriteDerivedOutput = 1; % OR inputs.OverwriteDerivedOutput = 1;
inputs.BiasCorr = 0; % OR inputs.BiasCorr = 0; 


%% Select dataset(s) to use
inputs.DataType = {'UKCP18'}; % {'UKCP18','UKCP18prob','HadUKGrid','ERA5','CMIP6'};
inputs.Variable = {'ATmin'}; % prefixes: 'T', 'VP', 'sWBGT', 'HD', 'AT', 'THI'; suffixes: 'max', 'min' or nothing for daily mean
inputs.Dataset = {'GCM-01','GCM-04','GCM-05'};
inputs.MMM = 1; % OR inputs.MMM = 1; % Calculate multi-model mean of all selected datasets
inputs.MMP = []; % OR inputs.MMP = 0.1 - 99.9; % Calculate multi-model percentile of all selected datasets


%% Subsetting of data for analysis
inputs.TemporalRange = [19810101, 20001230]; % yyyy, yyyymm or yyyymmdd start and end dates
inputs.AnnSummer = 'JJA'; % OR inputs.AnnSummer = 'Summer'; 'JJA'; 'MJJAS','ann'; % Summer = 1st June-15th Sept. or leave blank for annual mean
inputs.SpatialRange = [52,57;-5, 0]; % [latS,latN;lonW,lonE] for boxed region or [lat,lon] for single grid point
inputs.Region = []; % 'Scotland','North East','North West','Yorkshire and 
    % the Humber','East Midlands','West Midlands','East of England','Greater
    % London','South East','South West','Wales','Northern Ireland','Isle of
    % Man','RoI','Ireland','GB'

    
    
    
    
% inputs.TempRes = 'daily'; % OR inputs.TempRes = 'monthly';
% inputs.Resolution = {'2km'};
% inputs.CPM_period = 1981; % 1981 OR 2021 OR 2061



inputs.Years = [1981 2000];
inputs.AnnSummer = 'JJA'; % OR inputs1.AnnSummer = 'Summer'; 'JJA'; 'MJJAS'; % Summer = 1st June-15th Sept.
inputs.ExtremeMeanPctile = [95 99];
inputs.Pctile = [95 99];
% inputs2.DDx = 66;
inputs.OutputType = {'map'};
% inputs2.
inputs.Region = 'all'; % 'all' for grid cell output OR 'Scotland','North East','North West','Yorkshire and the Humber','East Midlands','West Midlands','East of England','Greater London','South East','South West','Wales','Northern Ireland','Isle of Man'


