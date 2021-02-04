function [data,xyz] = HEAT_step1(inputs,DataType,Dataset,Variable)
% [data,xyz] = HEAT_step1(inputs1,DataType,Dataset,Variable)
%
% Step 1 of the HEAT analysis code. This function loads the appropriate
% dataset and saves a derived variable if necessary.
%
% Outputs of data (the climate variable) and xyz (lat, long, time) can be
% produced by the function if the user wants to manually inspect the data
% that has been saved, however this is not necessary.


%% Set defaults
% Set directory paths
init_HEAT

% % For loading CPM data, a time period must be selected. This is 1981-2000
% % by default, but will not affect the loading of the RCM or GCM
% % slimulations
% if isfield(inputs,'CPM_period')
%     if inputs.CPM_period > 1979 && inputs.CPM_period < 2001
%         CPM_period = [1981,1991];
%     else
%         if inputs.CPM_period > 2019 && inputs.CPM_period < 2041
%             CPM_period = [2021,2031];
%         else
%             if inputs.CPM_period > 2059 && inputs.CPM_period < 2081
%                 CPM_period = [2061,2071];
%             end
%         end
%     end
%
% else
%     CPM_period = 1981;
% end


% Find if loading daily mean, min or max
if strcmp(Variable(end),'x')
    Tvar = 'tasmax';
elseif strcmp(Variable(end),'n')
    Tvar = 'tasmin';
else
    Tvar = 'tas';
end


% Set file name for derived variable based on inputs
fname = [Variable,'-',char(inputs.Domain),'-',Dataset];
disp('Checking if this data/variable combination has been loaded before:')


%% Check if this file has already been derived:
froot = [Deriveddir,fname]; % Take the file name...
files = dir([froot '*.nc']); % Then check if any files exist with this root

% If not, then load and process accordingly, then save
if isempty(files)
    disp(['No existing derived data file in ',Deriveddir])
    
    % If variable has been derived previously, see if overwriting is required
elseif inputs.OverwriteDerivedOutput == 1
    
    % Remove existing file in this case
    disp('Overwriting existing derived data file')
    for f = 1:length(files)
        file = [files(f).folder,'/',files(f).name];
        delete(file)
    end
else
    disp(['Derived data already exists in ',Deriveddir]')
end


%% Start loading data
disp('Processing derived data')

% Load elevation and lat-long data
load_xyz


if strcmp(DataType,'UKCP18')
    
    % Load  required simulation
    disp(['-> ',Dataset])
    
    % Change simulation name to UKCP18 file format and select elevation data at correct resolution
    if strcmp(Dataset(1:2),'RC')
        model = ['rcm85',Dataset(5:6)];
        ht = ht_12km;
        res = '12km/';
        template = [UKCP18dir,'12km/',Tvar,'/run',Dataset(5:6),'/',Tvar,'_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),'_day_19801201-19901130.nc'];
        
    elseif strcmp(Dataset(1:2),'CP')
        model = ['cpm85',Dataset(5:6)];
        ht = ht_2km;
        res = '2km/';
        template = [UKCP18dir,'2km/',Tvar,'/run',Dataset(5:6),'/',Tvar,'_rcp',model(4:5),'_land-cpm_uk_2.2km_',model(6:7),'_day_19801201-19811130.nc'];
        
    elseif strcmp(Dataset(1:2),'GC')
        model = ['gcm85',Dataset(5:6)];
        ht = ht_60km;
        res = '60km/';
        template = [UKCP18dir,'60km/',Tvar,'/run',Dataset(5:6),'/',Tvar,'_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),'_day_19791201-19891130.nc'];
        
    elseif strcmp(Dataset(1:2),'CM')
        model = ['gcm85',Dataset(7:8)];
        ht = ht_60km;
        res = '60km/';
        template = [UKCP18dir,'60km/',Tvar,'/run',Dataset(7:8),'/',Tvar,'_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),'_day_19791201-19891130.nc'];
    end
    
    % Directory of raw data for each required variable
    Tdir = [UKCP18dir,res,Tvar,'/run',Dataset(5:6),'/'];
    PSLdir = [UKCP18dir,res,'psl','/run',Dataset(5:6),'/'];
    SHdir = [UKCP18dir,res,'huss','/run',Dataset(5:6),'/'];
    
    % Find how many files are to be loaded/produced
    Tfiles = dir([Tdir '*.nc']);
    PSLfiles = dir([PSLdir '*.nc']);
    SHfiles = dir([SHdir '*.nc']);
    
    % Check if temporally consistent files are available for all variables
    if length(Tfiles) ~= length(PSLfiles) ~= length(SHfiles)
        
        % If not, find which time steps are
        matchingtimes = nan(length(Tfiles),3);
        
        % Go through each file for each variable to see which correspond to
        % each other
        for i = 1:length(Tfiles)
            for j = 1:length(PSLfiles)
                for k = 1:length(SHfiles)
                    if strcmp(Tfiles(i).name(end-19:end),PSLfiles(j).name(end-19:end)) && strcmp(Tfiles(i).name(end-19:end),SHfiles(k).name(end-19:end))
                        matchingtimes(i,1) = i;
                        matchingtimes(i,2) = j;
                        matchingtimes(i,3) = k;
                    end
                end
            end
        end
        
        % Keep only the time steps which are there for all variables
        matchingtimes = matchingtimes(~isnan(matchingtimes(:,1)),:);
        Tfiles = Tfiles(matchingtimes(:,1));
        PSLfiles = PSLfiles(matchingtimes(:,2));
        SHfiles = SHfiles(matchingtimes(:,3));
        
    end
    
    % Load required variable
    disp(['---> ',Variable])
    
    
    % Go through each file
    for f = 1:length(Tfiles)
        
        % Get the lists of files to load for each variable
        Tfile = [Tfiles(f).folder,'/',Tfiles(f).name];
        PSLfile = [PSLfiles(f).folder,'/',PSLfiles(f).name];
        SHfile = [SHfiles(f).folder,'/',SHfiles(f).name];
        
        % Load temperature
        T = ncread(Tfile,Tvar);
        ymd1 = ncread(Tfile,'yyyymmdd');
        time1 = ncread(Tfile,'time');
        
        % Load pressure and humidity to calculate vapour pressure
        PSL = ncread(PSLfile,'psl');
        % Convert pressure units in RCM simulations
        if strcmp(model(1),'r')
            PSL = PSL/100;
        end
        ymd2 = ncread(PSLfile,'yyyymmdd');
        [T,ymd,PSL,~] = check_consistent_timestep(T,ymd1,PSL,ymd2);
        P = p_surf(PSL,T,ht);
        SH = ncread(SHfile,'huss');
        ymd4 = ncread(SHfile,'yyyymmdd');
        [P,ymd,SH,~] = check_consistent_timestep(P,ymd,SH,ymd4);
        VP = VapourPressure(SH,P);
        
        % Calculate the correct heat stress metric
        if strcmp(Variable(1:2),'sW')
            data = SWBGTVP(T,VP);
            
        elseif strcmp(Variable(1:2),'HD')
            data = HumidexVP(T,VP);
            
        elseif strcmp(Variable(1:2),'AT')
            data = AppTempVP(T,VP);
            
        elseif strcmp(Variable(1:2),'VP')
            data = VP;
        end
        
        % Store date info
        [time1,~,~,~] = check_consistent_timestep(time1,ymd1,ymd,ymd);
        xyz.dates = ymd;
        xyz.time = time1;
        
        % Save the derived variable
        save_derived_nc(fname,data,xyz,Variable,template)
        
    end
end
