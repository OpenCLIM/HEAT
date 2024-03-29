% inputs_files_DAFNI.m
% 
% This .m file generates the template input file running HEAT on DAFNI. The
% majority of inputs are provided as environment variables.
% 

%% Set some basic experiment details
inputs.ExptName = 'HEAT_DAFNI_default';
inputs.OverwriteExpt = 1; % OR inputs.OverwriteExpt = 1;


%% Update fields if Environment variables are provided
disp('Updating inputs with environment variables')
% Then overwrite defaults with environment variables if running on DAFNI:
env_expn = getenv('EXPNAME');
env_runn = getenv('RUNNAME');
env_varn = getenv('VARNAME');
env_varno = getenv('VARNAMEOTHER');
env_scen = getenv('SCENARIO');
env_bases = getenv('BASELINE_S');
env_tims = getenv('TIMEPERIOD_S');
env_timl = getenv('TIMEPERIOD_L');
env_cal = getenv('CALENDAR');
env_summer = getenv('SUMMER');
env_outp = getenv('OUTPUT');
env_adap = getenv('ADAPT');
env_green = getenv('GREENING');
env_uhii = getenv('UHI_I');
env_exmt = getenv('EXMEAN');
env_abst = getenv('ABS_T');
env_pert = getenv('PER_T');
env_ddx = getenv('DDX');
env_x1 = getenv('X1');
env_x2 = getenv('X2');
env_y1 = getenv('Y1');
env_y2 = getenv('Y2');
env_reg = getenv('REGION');

% Experiment name for saving output
if ~isempty(env_expn)
    disp(['Environment variable found for Experiment Name: updating inputs file to ', char(string(env_expn))])
    inputs.ExptName = char(string(env_expn));
end

% Specify which UKCP18 member is being used to ensure correct time period
% is loaded corresponding to a warming level
if ~isempty(env_runn)
    disp(['Simulation name specified (for selecting correct time period): updating inputs file to ', char(string(env_runn))])
    inputs.Dataset = {env_runn};
end

% Select what climate variable to load and process
if ~isempty(env_varn)
    disp('Environment variable found for Variable: updating inputs file')
    inputs.Variable = {env_varn};
    % If user has selected 'Other', take the alternative name (for use for
    % example if temperature has a non standard name, e.g. Tmean)
    if strcmp(inputs.Variable,'Other')
        inputs.Variable = {env_varno};
    end
else % For offline testing
    inputs.Variable = {'tas'};
end
disp(['Using variable name ', char(string(inputs.Variable))])

% If user wants to select a warming level above pre-industrial
if ~isempty(env_scen)
    disp(['Environment variable found for Scenario: updating inputs file to ', char(string(env_scen))])
    inputs.Scenario = string(env_scen);
end

% Define a baseline start for acclimatisation
if ~isempty(env_bases)
    disp(['Environment variable found for defining baseline period: updating start year to ', char(string(env_bases))])
    inputs.BaselineStart = str2double(string(env_bases));
end

% Alternatively, user can specify a time period based on its start year...
if ~isempty(env_tims)
    disp(['Environment variable found for defining time period: updating start year to ', char(string(env_tims))])
    inputs.PeriodStart = str2double(string(env_tims));
end

% ...and the length of the period
if ~isempty(env_timl)
    inputs.PeriodLength = str2double(string(env_timl));
else
    inputs.PeriodLength = 30; % Otherwise assume 30 year default
end

% If user wants to define which calendar to use (360 day, 365 day, no leap)
if ~isempty(env_cal)
    disp('Environment variable found for setting calendar length: updating')
    inputs.Calendar = char(string(env_cal));
else 
    inputs.Calendar = 'd360'; % Default assumes 360 day calendar for UKCP18
end

% If user wants to use just summer period, not annual
if ~isempty(env_summer)
    disp('Environment variable found specifying analysis for summer/annual: updating')
    inputs.AnnSummer = char(string(env_summer));
else 
    inputs.AnnSummer = 'Ann'; % Default assumes 360 day calendar for UKCP18
end


% If user wants to use just summer period, not annual
if ~isempty(env_reg)
    disp('Environment variable found specifying UK region : updating')
    inputs.Region = char(string(env_reg));
end

% If a spatial subset has been specified, setup array to store info
if ~isempty(env_x1) && ~isempty(env_x2) % if range is given
        inputs.SpatialRange = nan(2,2);
elseif ~isempty(env_x1) % if single value is given
        inputs.SpatialRange = nan(2,1);
end

% Set spatial limits
if ~isempty(env_x1)
    disp('Environment variable found for defining x limit: updating inputs file')
    inputs.SpatialRange(1,1) = str2double(string(env_x1));
end
if ~isempty(env_x2)
    disp('Environment variable found for defining second x limit: updating inputs file')
    inputs.SpatialRange(1,2) = str2double(string(env_x2));
end
if ~isempty(env_y1)
    disp('Environment variable found for defining y limit: updating inputs file')
    inputs.SpatialRange(2,1) = str2double(string(env_y1));
end
if ~isempty(env_y2)
    disp('Environment variable found for defining second y limit: updating inputs file')
    inputs.SpatialRange(2,2) = str2double(string(env_y2));
end

% Choose the output type (e.g. a netCDF for use in HARM, some analysis etc.)
if ~isempty(env_outp)
    disp(['Environment variable found for output type: carrying out analysis for  to ', char(string(env_outp))])
    inputs.OutputType = env_outp;
end

% Choose a percentile for simulating best-case-scenario acclimatisation ?
% HARM MMT and RR will be adjusted in line with this warming
if ~isempty(env_adap)
    disp('Environment variable found for setting acclimatisation percentile: updating inputs file')
    inputs.MMTpctile = str2double(string(env_adap));
end


% Set the urban greening/cooling potential (how many degrees celsius fully
% urbanised grid cells will have local temperatures reduced by)
if ~isempty(env_green)
    disp('Environment variable found for parameterising cooling due to greening: updating inputs file')
    inputs.Greening = str2double(string(env_green));
else
    inputs.Greening = 0;
end

% Set the Urban Heat Island Intensity (how many degrees celsius full
% urbanisation will increase local temperatures)
if ~isempty(env_uhii)
    disp('Environment variable found for defining Urban Heat Island intensity: updating inputs file')
    inputs.UHI_I = str2double(string(env_uhii));
else
    inputs.UHI_I = 2;
end

% Set extreme mean threshold
if ~isempty(env_exmt)
    disp('Environment variable found for setting extreme mean percentile threshold: updating inputs file')
    inputs.ExtremeMeanPctile = str2double(string(env_exmt));
end

% Set an absolute threshold for extremes analysis
if ~isempty(env_abst)
    disp('Environment variable found for setting absolute threshold: updating inputs file')
    inputs.AbsThresh = str2double(string(env_abst));
end

% Choose a percentile threshold for extremes analysis
if ~isempty(env_pert)
    disp('Environment variable found for setting percentile threshold: updating inputs file')
    inputs.PercentileThresh = str2double(string(env_pert));
end

% Choose a percentile above which to calculate degree day metric
if ~isempty(env_ddx)
    disp('Environment variable found for setting degree day percentile: updating inputs file')
    inputs.DD = str2double(string(env_ddx));
end


disp(' ')

