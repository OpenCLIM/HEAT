function [] = HEAT(inputs,varargin)
% HEAT v.1.0
%
% Run the OpenCLIM Heat Extremes Analysis Toolbox (HEAT). This function
% does an inital setup and check of file directories, then goes through
% each required dataset, loading, processing and generating output
% accordingly.
%
% Input files, in the form of structures, can be created using the format
% provided in launch_file_TEMPLATE.m or interactively using
% generate_launch_file.m. Alternatively the inputs can be defined directly
% from a script (e.g. launch_file_TEMPLATE.m) if running HEAT app
% standalone.
%
% inputs = information regarding data which is to be loaded, including
% variables, dataset to use, temporal and spatial domain, experiment name
% and whether to save derived variables.
%
% Coded by A.T. Kennedy-Asser, University of Bristol, 2021.
% Contact: alan.kennedy@bristol.ac.uk
%

%% Initialise
disp('Running HEAT v.1.0')
disp('-----')


% Set directory paths
init_HEAT

% Record start time
startt = now;


%% Check inputs
% If running HEAT as standalone app, inputs should be included in a
% seperate inputs script, which can be adapted from input_files_TEMPLATE.m.
% If running HEAT through the MATLAB GUI, thne it is possible to have run
% input_files_TEMPLATE.m first, in which case the structures this script
% creates can be input.

% Check if inputs is a script, in which case run it
if ~exist('inputs','var')
    run('input_files_DAFNI.m')
end

if ischar(inputs)
    if exist(inputs,'file')
        run(inputs)
        
    else
        disp('Error: input script not found');
        return
    end
    
    % Otherwise if inputs is a structure, then adjust the name ready to run HEAT
else
    if ~isstruct(inputs)
        disp('Error: inputs1 must be a structure');
        return
    end
    if ~isfield(inputs,'ExptName')
        disp('Error: inputs appears to be the wrong structure');
        return
    else
        inputs = inputs;
    end
end

% Set default domain to UK if not specified
if ~isfield(inputs,'Domain')
    inputs.Domain = 'UK';
end

%% Find which steps of HEAT to run
% Set default to not run steps
runstep1 = 0;
runstep2 = 0;
runstep3 = 0;

% Run steps if necessary
if isfield(inputs,'SaveDerivedOutput')
    if inputs.SaveDerivedOutput == 1
        runstep1 = 1;
    end
end

if isfield(inputs,'OutputType')
    runstep2 = 1;
end

if isfield(inputs,'WorkflowOutput')
    runstep3 = 1;
end



% % If multiple inputs are provided, go through to find if running further steps:
% for i = 1:length(varargin)
%
%     % Set inputs2 and inputs3 if they are provided
%     if isstruct(varargin{i})
%         if isfield(varargin{i},'OutputType')
%             inputs2 = varargin{i};
%             % Switch on step
%             runstep2 = 1;
%         else
%             if isfield(varargin{i},'ExampleAdaptationParameter')
%                 inputs3 = varargin{i};
%                 % Switch on step
%                 runstep3 = 1;
%             end
%
%             % Can add further steps here in the future as necessary.
%         end
%     end
% end

% % If inputs were provided as a script, structures inputs2, inputs3 etc. may
% % exist, in which case these steps need run:
% if exist('inputs2','var')
%     if isstruct(inputs2)
%         runstep2 = 1;
%     end
% end
%
% if exist('inputs3','var')
%     if isstruct(inputs3)
%         runstep3 = 1;
%     end
% end


%% Set up output directory
% First check the experiment name won't overwrite the DerivedData directory
if strcmp(inputs.ExptName,'DerivedData')
    disp('Cannot call experiment "DerivedData": CANCELLING')
    return
end

% Next, check if output directory already exists
if exist([Outputdir,'/',inputs.ExptName],'dir')
    
    % Check if overwriting has been authorised in input file
    if ~isfield(inputs,'OverwriteExpt')
        disp('Permission to overwrite unclear in input file: CANCELLING')
        disp('-----')
        % Cancel the run
        return
    else
        
        % If it does, check whether or not to overwrite
        disp('Warning: existing experiment exists with this name')
        if inputs.OverwriteExpt == 1
            disp('Overwriting enabled in input file: Re-running experiment')
            disp('-----')
            
        else
            disp('Overwriting disabled in input file: CANCELLING')
            disp('-----')
            % If not overwriting, cancel the run
            return
        end
        
    end
    
else
    % Create output directory
    mkdir([Outputdir,'/',inputs.ExptName])
end


% Save input files for future reference
save([Outputdir,'/',inputs.ExptName,'/inputs.mat'],'inputs')


%% Produce derived data if required (step 1)
if runstep1 == 1
    
    % Go through each required dataset/simulation/variable
    % Data hierarchy:
    % DataType (e.g. UKCP18, ERA5, CMIP6) -> Dataset (e.g. specific simulation, observation resolution)
    
    for d = 1:length(inputs.DataType)
        DataType = char(inputs.DataType(d));
        
        % Load model data if required
        if ismember(inputs.DataType(d),'UKCP18')
            
            % Load each required simulation
            for s = 1:length(inputs.Dataset)
                Dataset = char(inputs.Dataset(s));
                
                % Only attempt to load UKCP18-related data
                if strcmp(Dataset(1),'G') || strcmp(Dataset(1),'R') || strcmp(Dataset(1),'C')
                    
                    % Load each required variable
                    for v = 1:length(inputs.Variable)
                        Variable = char(inputs.Variable(v));
                        
                        if strcmp(Variable(1),'T') % Any temperature variable might as well be loaded from raw data for UKCP18 models
                            disp('No benefit in saving derived data: simply use raw temperature data')
                            disp('-----')
                        else
                            % Run Step 1: Produce derived data
                            HEAT_step1(inputs,DataType,Dataset,Variable);
                        end
                    end
                end
            end
            
        end
        
        
        % Load observational data if required
        if ismember(inputs.DataType(d),'HadUKGrid')
            
            % Load each required simulation
            for s = 1:length(inputs.Dataset)
                Dataset = char(inputs.Dataset(s));
                
                % Only attempt to load UKCP18-related data
                if strcmp(Dataset(1),'1') || strcmp(Dataset(1),'2') || strcmp(Dataset(1),'6')
                    
                    % Load each required variable
                    for v = 1:length(inputs.Variable)
                        Variable = char(inputs.Variable(v));
                        
                        % If using 1km, 12km or 60km there is no benefit in
                        % deriving anything for Tmax or Tmin
                        if ~strcmp(Dataset,'2km')
                            if strcmp(Variable,'Tmin') || strcmp(Variable,'Tmax') % Tmax and Tmin might as well be loaded from raw data for UKCP18 models, however Tmean could be derived at daily resolution
                                disp('No benefit in saving derived data: simply use raw max./min. temperature data')
                                disp('-----')
                            elseif strcmp(Variable,'VPmax') || strcmp(Variable,'VPmin')
                                disp('No benefit in saving derived data: all VP data is daily mean')
                                disp('-----')
                            else
                                % Run Step 1: Produce derived data
                                HEAT_step1(inputs,DataType,Dataset,Variable);
                            end
                            
                            % HadUK-Grid is not available at 2km resolution so can
                            % be regridded as derived data for any variable
                        else
                            % Run Step 1: Produce derived data
                            HEAT_step1(inputs,DataType,Dataset,Variable);
                        end
                        
                    end
                end
            end
            
        end
        
    end
end


%% Produce diagnostic or workflow data if required (steps 2 and 3)
if runstep2 == 1
    
    % Go through each required dataset/simulation/variable
    % Data hierarchy:
    % DataType (e.g. UKCP18, ERA5, CMIP6) -> Dataset (e.g. specific simulation, observation resolution)
    
    % Load each required variable
    for v = 1:length(inputs.Variable)
        Variable = char(inputs.Variable(v));
        
        % Run Step 2: Extremes analysis
        if runstep2 == 1
            HEAT_step2(inputs,Variable)
        end
        
    end
  
        
%         % Load model data if required
%         if ismember(inputs.DataType(d),'HadUKGrid')
%             
%             % Load each required simulation
%             for s = 1:length(inputs.Dataset)
%                 Dataset = char(inputs.Dataset(s));
%                 
%                 % Load each required variable
%                 for v = 1:length(inputs.Variable)
%                     Variable = char(inputs.Variable(v));
%                     
%   
%                     %% Run Step 2: Extremes analysis
%                     if runstep2 == 1
%                         HEAT_step2(inputs,data,xyz,Dataset,Variable,inputs.ExptName)
%                     end
%                     
%                 end
%             end
%             
%         end
        
%     end
end


%% Generate output for other models in workflows
if runstep3 == 1
    
    % Check if this file has already been derived:
    froot = [Outputdir,'/',inputs.ExptName,'/']; % Take the file name...
    files = dir([froot '*',Variable,'.csv']); % Then check if any files exist with this root
    
    % If so, delete
    if ~isempty(files)
        for f = 1:length(files)
            file = [files(f).folder,'/',files(f).name];
            delete(file)
        end
    end
    
    for d = 1:length(inputs.DataType)
        DataType = char(inputs.DataType(d));
        
        % Load each required variable
        for v = 1:length(inputs.Variable)
            Variable = char(inputs.Variable(v));
            
            HEAT_step3(inputs,DataType,Variable);
            
        end
    end
    
end


%% Finish up
disp(' ')
disp(['HEAT run "',inputs.ExptName,'" complete',])
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))
disp('-----')


