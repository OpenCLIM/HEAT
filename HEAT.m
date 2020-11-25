function [] = HEAT(inputs1,varargin)
% HEAT v.1.0
% 
% Run the OpenCLIM Heat Extremes Analysis Toolbox (HEAT) for selected input
% data. 
% 
% Input files, in the form of structures, can be created using the format
% provided in launch_file_TEMPLATE.m or interactively using
% generate_launch_file.m.
% 
% inputs1 = information regarding data which is to be loaded, including
% variables, dataset to use, temporal and spatial domain, experiment name
% and whether to save derived variables.
% 
% inputs2 = specific output for analysis, including thresholds for extremes
% analysis, plotting information.
% 
% inputs3,...,x = options for additional modules accounting for
% socio-economic data, adaptation, agriculture etc. (NOT YET INCLUDED).
% 
% Coded by A.T. Kennedy-Asser, University of Bristol, 2020.
% 

%% Initialise %%
disp('Running HEAT v.1.0')
disp('-----')

% Set directory paths
init_HEAT

% Record start time
startt = now;


%% Check inputs %%
% Check inputs1 is correct (it is essential to run)
if ~isstruct(inputs1)
    disp('Error: inputs1 must be a structure');
    return
end
if ~isfield(inputs1,'Domain')
    disp('Error: inputs1 appears to be the wrong structure');
    return
end

% Go through other input variables to find if running further steps:
% Set switch to not run additional steps as default
csvout = 0;
runstep2 = 0;
runstep3 = 0;

for i = 1:length(varargin)
    
    % Set inputs2 and inputs3 if they are provided
    if isstruct(varargin{i})
        if isfield(varargin{i},'OutputType')
            inputs2 = varargin{i};
            % Switch on step
            runstep2 = 1;
        else
            if isfield(varargin{i},'ExampleAdaptationParameter')
                inputs3 = varargin{i};
                % Switch on step
                runstep3 = 1;
            end
        end
    end
end


%% Set up output directory %%
% First, check if output directory already exists
if exist([Deriveddir,'/',inputs1.ExptName],'dir')
    
    % If it does, check whether or not to overwrite
    disp('Warning: existing experiment exists with this name.')
    if inputs1.OverwriteDerivedOutput == 1
        disp('Overwriting enabled in input file')
        disp('-----')

    else
        disp('Overwriting disabled in input file: CANCELLING')
        disp('-----')
        % If not overwriting, cancel the run
        return
    end
    
else % Create output directory
    mkdir([Deriveddir,'/',inputs1.ExptName])
end

% Save input files for future reference
save([Deriveddir,'/',inputs1.ExptName,'/inputs1.mat'],'inputs1')
if exist('inputs2','var')
    save([Deriveddir,'/',inputs1.ExptName,'/inputs2.mat'],'inputs2')
end
if exist('inputs3','var')
    save([Deriveddir,'/',inputs1.ExptName,'/inputs3.mat'],'inputs3')
end


%% Load each required dataset/simulation/variable %%
for d = 1:length(inputs1.DataType)
    DataType = char(inputs1.DataType(d));
    
    % Load model data if required
    if ismember(inputs1.DataType(d),'model')
        
        % Load each required simulation
        for s = 1:length(inputs1.Simulation)
            Simulation = char(inputs1.Simulation(s));
            
            % Load each required variable
            for v = 1:length(inputs1.Variable)
                Variable = char(inputs1.Variable(v));
                
                % Set file name for derived variable based on inputs
                fname = ['HEAT_v1.0-',char(inputs1.Domain),'-',Simulation,'-',Variable,'.nc'];
                disp('Checking if this data/variable combination has been loaded before:')
                
                % Check if this file has already been derived:
                % If not, then load and process accordingly, then save
                if ~exist([Deriveddir,fname],'file')
                    disp(['No existing derived data file in ',Deriveddir])
                    
                    % Load the raw data
                    [data,xyz,template] = load_data(DataType,Simulation,Variable);
                    
                    % Save the derived variable if requested
                    % (Note, Tmean, Tmin and Tmax are not saved as the derived variable is no quicker to load than the raw data)
                    if inputs1.SaveDerivedOutput == 1
                        if ~strcmp(Variable(1:2),'Tm')
                            save_derived_nc(fname,data,xyz,Variable,template)
                        end
                    end
                    
                    
                else % If variable has been derived previously, see if overwriting is required
                    if inputs1.OverwriteDerivedOutput == 1
                        
                        % Remove existing file in this case
                        disp('Overwriting existing derived data file')
                        delete([Deriveddir,fname])
                        
                        % Load the raw data
                        [data,xyz,template] = load_data(DataType,Simulation,Variable);
                        
                        % Save the derived variable if requested
                        if inputs1.SaveDerivedOutput == 1
                            if ~strcmp(Variable(1:2),'Tm')
                                save_derived_nc(fname,data,xyz,Variable,template)
                            end
                        end

                        
                    else % If the derived data exists and it is okay to use, load it
                        disp('Derived variable has already been calculated: Loading')
                        data = ncread([Deriveddir,fname],Variable);
                        xyz.dates = ncread([Deriveddir,fname],'yyyymmdd');
                        

                    end
                end
                
                
                
                %% Run Step 2: Extremes analysis %%
                if runstep2 == 1
                    disp('Running extremes analysis')
                    
                    % Temporally subset data as required
                    [data,xyz.dates] = subset_temporal(data,xyz.dates,inputs2.Years,inputs2.AnnSummer);
                    
                    % Spatially subset data as required
                    
                    
                    % Calculate percentile if required
                    if isfield(inputs2,'Pctile')
                        % Create empty array for output
                        data_pctiles = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs2.Pctile));
                        
                        for p = 1:length(inputs2.Pctile)
                            data_pctiles(:,:,p) = prctile(data,inputs2.Pctile(p),3);
                        end
                    end
                    
                    % Calculate extreme mean if required
                    if isfield(inputs2,'ExtremeMeanPctile')
                        % Create empty array for output
                        data_ExMean = nan(length(data(:,1,1)),length(data(1,:,1)),length(inputs2.ExtremeMeanPctile));
                        
                        for p = 1:length(inputs2.ExtremeMeanPctile)
                            
                            Txx_temp = data >= prctile(data,inputs2.ExtremeMeanPctile(p),3);
                            Txx = nan(size(Txx_temp));
                            Txx(Txx_temp == 1) = 1;
                            data_Txx = data .* Txx;
                            data_ExMean(:,:,p) = squeeze(nanmean(data_Txx,3));
                        end
                    end
                    
                    
                    % Generate requested output
                    disp('Producing output')
                    for o = 1:length(inputs2.OutputType)
                        OutputType = string(inputs2.OutputType(o));
                        
                        % Generate map if requested
                        if strcmp(OutputType,'map')
                            figure
                            UK_subplot(data_ExMean,[Simulation,' ',Variable, ' extreme mean >95th percentile'])
%                             figure
%                             UK_subplot(data_pctiles,[Simulation,' ',Variable, ' 95th percentile'])
                        end
                        
                        
                        % Generate ARCADIA input csv if requested
%                         csvout = 0;
                        if strcmp(OutputType,'ARCADIA')
                            csvout = 1;
                            
                            if ~exist('data_csv','var')
                                data_csv.(Variable) = data;
                            else
                                if isfield(data_csv,Variable)
                                    
                                    data_csv.(Variable) = cat(4,data_csv.(Variable),data);
                                else
                                    data_csv.(Variable) = data;
                                end
                            end
                            
                            
                            
                        end
                        
                    end
                    disp('-----')
                    
                end
                
                %% Run Step 3: Produce output %%
                if runstep3 == 1
                    
                    
                    
                end
                
                
                
                
            end
        end
        
        % Save csv if required
        if csvout == 1
            
            for v = 1:length(inputs1.Variable)
                Variable = char(inputs1.Variable(v));
                
                for i = 1:length(data_csv.(Variable)(:,1,1,1))
                    for j = 1:length(data_csv.(Variable)(1,:,1,1))
                        %                 for k = 1:length(data_csv.sWBGT(1,1,1,:))
                        
                        data_csv_temp = squeeze(data_csv.(Variable)(i,j,:,:))';
                        
                        csv_name = [Deriveddir,'/',inputs1.ExptName,'/',num2str(i),num2str(j),'_test_',Variable,'.csv'];
                        %                     writematrix(data_csv_temp,csv_name)
                        csvwrite(csv_name,data_csv_temp)
                        %                 end
                    end
                end
            end
        end
    end
end

disp(['HEAT run "',inputs1.ExptName,'" complete',])
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))
disp('-----')


