function [] = HEAT(inputs,varargin)
% HEAT v.1.0
% 
% Run the OpenCLIM Heat Extremes Analysis Toolbox (HEAT). This function
% does an inital setup and check of file directories, then goes through
% each required dataset, loading (step 1), processing and generating output
% (step 2) accordingly.
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
% inputs2 = specific output for analysis, including thresholds for extremes
% analysis, plotting information.
% 
% inputs3,...,x = options for additional modules accounting for
% socio-economic data, adaptation, agriculture etc. (NOT YET INCLUDED).
% 
% Coded by A.T. Kennedy-Asser, University of Bristol, 2020.
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
if isstr(inputs)
    if exist(inputs) == 2
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
    if ~isfield(inputs,'Domain')
        disp('Error: inputs1 appears to be the wrong structure');
        return
    else
        inputs1 = inputs;
    end
end

%% Find which steps of HEAT to run
% Note: Step 1 is essential and is always run.

% Set default to not run extra steps
runstep2 = 0;
runstep3 = 0;

% If multiple inputs are provided, go through to find if running further steps:
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
            
            % Can add further steps here in the future as necessary.
        end
    end
end

% If inputs were provided as a script, structures inputs2, inputs3 etc. may
% exist, in which case these steps need run:
if exist('inputs2','var')
    if isstruct(inputs2)
        runstep2 = 1;
    end
end

if exist('inputs3','var')
    if isstruct(inputs3)
        runstep3 = 1;
    end
end


%% Set up output directory 
% First check the experiment name won't overwrite the DerivedData directory
if strcmp(inputs1.ExptName,'DerivedData')
    disp('Cannot call experiment "DerivedData": CANCELLING')
    return
end

% Next, check if output directory already exists
if exist([Outputdir,'/',inputs1.ExptName],'dir')
    
    % If it does, check whether or not to overwrite
    disp('Warning: existing experiment exists with this name.')
    if inputs1.OverwriteExpt == 1
        disp('Overwriting enabled in input file')
        disp('-----')

    else
        disp('Overwriting disabled in input file: CANCELLING')
        disp('-----')
        % If not overwriting, cancel the run
        return
    end
    
else % Create output directory
    mkdir([Outputdir,'/',inputs1.ExptName])
end

% Save input files for future reference
save([Outputdir,'/',inputs1.ExptName,'/inputs1.mat'],'inputs1')
if exist('inputs2','var')
    save([Outputdir,'/',inputs1.ExptName,'/inputs2.mat'],'inputs2')
end
if exist('inputs3','var')
    save([Outputdir,'/',inputs1.ExptName,'/inputs3.mat'],'inputs3')
end


%% Go through each required dataset/simulation/variable
for d = 1:length(inputs1.DataType)
    DataType = char(inputs1.DataType(d));
    
    % Load model data if required
    if ismember(inputs1.DataType(d),'model')
        
        % Load each required simulation
        for s = 1:length(inputs1.Dataset)
            Dataset = char(inputs1.Dataset(s));
            
            % Load each required variable
            for v = 1:length(inputs1.Variable)
                Variable = char(inputs1.Variable(v));
                
                %% Run Step 1: Load data
                [data,xyz] = HEAT_step1(inputs1,DataType,Dataset,Variable);
                
                
                %% Run Step 2: Extremes analysis
                if runstep2 == 1
                    HEAT_step2(inputs2,data,xyz,Dataset,Variable,inputs1.ExptName)
                end
                
                %% Run Step 3: e.g. Adaptation modelling?
                if runstep3 == 1
                    
                end
                
                % Further steps can be added as the toolbox is developed
                
                
            end
        end
        
%         % Save csv if required after all models loaded
%         if csvout == 1
%             
%             for v = 1:length(inputs1.Variable)
%                 Variable = char(inputs1.Variable(v));
%                 
%                 for i = 1:length(data_csv.(Variable)(:,1,1,1))
%                     for j = 1:length(data_csv.(Variable)(1,:,1,1))
%                         %                 for k = 1:length(data_csv.sWBGT(1,1,1,:))
%                         
%                         data_csv_temp = squeeze(data_csv.(Variable)(i,j,:,:))';
%                         
%                         csv_name = [Outputdir,'/',inputs1.ExptName,'/',num2str(i),num2str(j),'_test_',Variable,'.csv'];
%                         %                     writematrix(data_csv_temp,csv_name)
%                         csvwrite(csv_name,data_csv_temp)
%                         %                 end
%                     end
%                 end
%             end
%         end
    end
end

disp(['HEAT run "',inputs1.ExptName,'" complete',])
endt = now;
fprintf('Total time taken to run: %s\n', datestr(endt-startt,'HH:MM:SS'))
disp('-----')


