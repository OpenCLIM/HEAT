function [data,xyz] = HEAT_step1(inputs1,DataType,Dataset,Variable)
% [data,xyz] = HEAT_step1(inputs1,DataType,Dataset,Variable)
% 
% Step 1 of the HEAT analysis code. This function loads the appropriate
% dataset and saves a derived variable if necessary.

%% Set defaults
% Set directory paths
init_HEAT

% For loading CPM data, a time period must be selected. This is 1981-2000
% by default, but will not affect the loading of the RCM or GCM
% slimulations
if isfield(inputs1,'CPM_period')
    CPM_period = inputs1.CPM_period;
else
    CPM_period = 1981;
end

% Set file name for derived variable based on inputs
fname = [Variable,'-',char(inputs1.Domain),'-',Dataset];
disp('Checking if this data/variable combination has been loaded before:')

% Check if this file has already been derived:
froot = [Deriveddir,fname]; % Take the file name...
file = dir([froot '*.nc']); % Then check if any files exist with this root

% If not, then load and process accordingly, then save
if isempty(file)
    disp(['No existing derived data file in ',Deriveddir])
    
    % Load the raw data
    [data,xyz,template] = load_data(DataType,Dataset,Variable,CPM_period);
    
    % Save the derived variable if requested
    % (Note, Tmean, Tmin and Tmax are not saved as the derived variable is no quicker to load than the raw data)
    if inputs1.SaveDerivedOutput == 1
        if ~strcmp(Variable(1:2),'Tm')
            save_derived_nc(fname,data,xyz,Variable,template)
        end
    end
    
    
else % If variable has been derived previously, see if overwriting is required
    if inputs1.SaveDerivedOutput == 1 && inputs1.OverwriteDerivedOutput == 1
        
        % Remove existing file in this case
        disp('Overwriting existing derived data file')
        file = [file.folder,'/',file.name];
        delete(file)
        
        % Load the raw data
        [data,xyz,template] = load_data(DataType,Dataset,Variable,CPM_period);
        
        % Save the derived variable if requested
        if inputs1.SaveDerivedOutput == 1
            if ~strcmp(Variable(1:2),'Tm')
                save_derived_nc(fname,data,xyz,Variable,template)
            end
        end
        
        
    else % If the derived data exists and it is okay to use, load it
        disp('Derived variable has already been calculated: Loading')
        disp(['-> Loading ',Dataset])
        disp(['--->  ',Variable])
        
        file = [file.folder,'/',file.name];
        data = ncread(file,Variable);
        xyz.dates = ncread(file,'yyyymmdd');
        disp('-----')
        
        
    end
end
