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

disp(' ')
disp('Running Step 1: producing derived data')
disp('-----')

% Find if loading daily mean, min or max
if strcmp(Variable(end),'x')
    Tvar = 'tasmax';
elseif strcmp(Variable(end),'n')
    Tvar = 'tasmin';
else
    Tvar = 'tas';
end

% Change the file name slightly for HadUK-Grid data
if strcmp(DataType,'HadUKGrid')
    Dataset1 = ['HadUK-Grid-',Dataset];
else
    Dataset1 = Dataset;
end

% Set file name for derived variable based on inputs
fname = [Variable,'-',char(inputs.Domain),'-',Dataset1];


%% Check if this file has already been derived:
disp('Checking if this data/variable combination has been loaded before:')
froot = [Deriveddir,fname]; % Take the file name...
files = dir([froot '*.nc']); % Then check if any files exist with this root

% If not, then load and process accordingly, then save
if isempty(files)
    disp(['No existing derived data file in ',Deriveddir])
    
    % If permission to overwrite derived data is not clearly requested, skip it
elseif ~isfield(inputs,'OverwriteDerivedOutput')
        disp('Permission to overwrite derived data unclear in input file: skipping')
        disp('-----')
        % Bybass the derived data processing
        skipload = 1;

    % If variable has been derived previously, see if overwriting is required
elseif inputs.OverwriteDerivedOutput == 1
    
    % Remove existing file in this case
    disp('Overwriting existing derived data file')
    for f = 1:length(files)
        file = [files(f).folder,'/',files(f).name];
        delete(file)
    end
else
    if ~exist('skipload','var')
        disp(['Derived data already exists in ',Deriveddir,': skipping'])
        disp('-----')
        skipload = 1;
    end
end


%% Start loading data
if ~exist('skipload','var')
    disp('Processing derived data')
    
    % Load elevation and lat-long data
    load_xyz
    
    %% UKCP18 data
    if strcmp(DataType,'UKCP18')
        
        % Load required simulation
        disp(['-> UKCP18 ',Dataset])
        
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
        
        % Check each variable for missing data files ? this was an issue
        % with some early UKCP18 products, although may be corrected in
        % future:
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
        
        %% HadUK-Grid data
        % The process for deriving HadUK-Grid data is slightly longer as it
        % also provides a method for deriving a daily mean temperature
        % netCDF file. This is the mean of the daily max. and daily min.
        % data, with the monthly mean of the daily values adjusted to that
        % it is consistent with the HadUK-Grid monthly daily mean
        % temperature data provided by the UK Met Office. This derived
        % daily mean temperature needs to be calculated prior to
        % calculating daily means of any of the heat stress variables.
        
    elseif strcmp(DataType,'HadUKGrid')
        
        % Load required resolution
        disp(['-> HadUK-Grid at ',Dataset])
        
        % Change simulation name to UKCP18 file format and select elevation data at correct resolution
        if strcmp(Dataset(1:2),'12')
            res = '12km';
            template = [HadUKdir,'tasmax/',res,'/','tasmax_hadukgrid_uk_',res,'_day_19600501-19600531.nc'];
            interpolate = 0;
        elseif strcmp(Dataset(1:2),'1k')
            res = '1km';
            template = [HadUKdir,'tasmax/1km/','tasmax_hadukgrid_uk_1km_day_19600501-19600531.nc'];
            interpolate = 0;
        elseif strcmp(Dataset(1:2),'60')
            res = '60km';
            template = [HadUKdir,'tasmax/',res,'/','tasmax_hadukgrid_uk_',res,'_day_19600501-19600531.nc'];
            interpolate = 0;
        elseif strcmp(Dataset(1:2),'2k')
            res = '1km';
            res2 = '2km';
            template = [HadUKdir,'tasmax/',res,'/','tasmax_hadukgrid_uk_',res,'_day_19600501-19600531.nc'];
            interpolate = 1;
        end
        
        % Directory of raw data for each required variable
        Tmaxdir = [HadUKdir,'tasmax/',res,'/'];
        Tmindir = [HadUKdir,'tasmin/',res,'/'];
        Tdir = [HadUKdir,'tas/',res,'/'];
        VPdir = [HadUKdir,'pv/',res,'/'];
        
        % Find how many files are to be loaded/produced
        Tmaxfiles = dir([Tmaxdir '*.nc']);
        Tminfiles = dir([Tmindir '*.nc']);
        
        % Check if temporally consistent files are available for all variables
        if length(Tmaxfiles) ~= length(Tminfiles)
%             disp('Different number of HadUK-Grid Tmax and Tmin files')
            
            % If not, find which time steps are
            matchingtimes = nan(length(Tmaxfiles),2);
            
            % Go through each file for each variable to see which correspond to
            % each other
            for i = 1:length(Tmaxfiles)
                for j = 1:length(Tminfiles)
                    if strcmp(Tmaxfiles(i).name(end-19:end),Tminfiles(j).name(end-19:end))
                        matchingtimes(i,1) = i;
                        matchingtimes(i,2) = j;
                    end
                end
            end
            
            % Keep only the time steps which are there for all variables
            matchingtimes = matchingtimes(~isnan(matchingtimes(:,1)),:);
            Tmaxfiles = Tmaxfiles(matchingtimes(:,1));
            Tminfiles = Tminfiles(matchingtimes(:,2));
        end
        
        
        % Find if deriving daily mean, min or max
        if strcmp(Variable(end),'x')
            Tvar = 'tasmax';
        elseif strcmp(Variable(end),'n')
            Tvar = 'tasmin';
        else
            Tvar = 'tas';
        end
        
        % HadUK-Grid VP cannot be loaded for daily max./min.
        if strcmp(Variable(1),'V')
            Tvar = 'tas';
        end
        
        disp(['---> ',Variable])
        
        % If deriving Tmean, use average of daily min. and max., then
        % adjust monthly mean to match monthly tas data
        if strcmp(Tvar,'tas')
            
            % Go through each file
            for f = 1:length(Tmaxfiles)
                
                % Get the lists of files to load for each variable
                Tmaxfile = [Tmaxfiles(f).folder,'/',Tmaxfiles(f).name];
                Tminfile = [Tminfiles(f).folder,'/',Tminfiles(f).name];
                
                % Find month and year that has been loaded for daily data
                month = str2double(Tmaxfile(end-6:end-5));
                year = Tmaxfile(end-10:end-7);
                
                % 1960 is not available for VP, so skip and start at 1961
                if str2double(year) > 1960
                    
                    % Check if Tmean has already been calculated
                    fname2 = ['T-',char(inputs.Domain),'-',Dataset1];
                    froot = [Deriveddir,fname2]; % Take the file name...
                    files = dir([froot,'-',Tmaxfile(end-19:end)]); % Then check if any files exist with this root
                    
                    % If not, then load and process accordingly, then save
                    if isempty(files) && ~strcmp(Variable(1),'V')
                        disp(['No existing Tmean derived data file in ',Deriveddir])
                        
                        % Load daily max and min temperatures
                        Tmax = ncread(Tmaxfile,'tasmax');
                        time1 = ncread(Tmaxfile,'time');
                        ymd1 = char(datetime(time1/24 + datenum(1800,01,01),'ConvertFrom','datenum','Format','yyyyMMdd'))';
                        Tmin = ncread(Tminfile,'tasmin');
                        
                        % Equivalent Tmean and vapour pressure file to load
                        Tfile = [Tdir,'tas_hadukgrid_uk_',res,'_mon_',year,'01-',year,'12.nc'];
                        Tmon = ncread(Tfile,'tas');
                        
                        % Subset the correct month
                        Tmon = Tmon(:,:,month);
                        
                        % Calculate Tmean as average of Tmax and Tmin
                        T = (Tmax + Tmin)/2;
                        
                        % Adjust monthly mean to equal the monthly data
                        T = T - (nanmean(T,3)-Tmon);
                        
                        % Store date info
                        xyz.dates = ymd1;
                        xyz.time = time1;
                        
                        % Convert from 1km to 2km resolution if required
                        if interpolate == 1
%                             startit = now;
                            T = change_res(T,res2);
%                             endit = now;
%                             fprintf('Total time taken to interpolate: %s\n', datestr(endit-startit,'HH:MM:SS'))
                        end
                        
                        % Save the derived variable
                        save_derived_nc(fname2,T,xyz,'T',template)
                        
                    end
                    
                    % Otherwise, if producing a heat stress variable or VP
                    if ~strcmp(Variable,'T')
                        
                        % Load Tmean if required (i.e. if derived data calculated by a previous run)
                        if ~strcmp(Variable(1:2),'VP') &&  ~exist('T','var')
                            T = ncread([Deriveddir,fname2,'-',Tmaxfile(end-19:end)],'T');
                        end
                        
                        % Load the appropriate VP file
                        VPfile = [VPdir,'pv_hadukgrid_uk_',res,'_mon_',year,'01-',year,'12.nc'];
                        VP = ncread(VPfile,'pv');
                        
                        % Subset the correct month
                        VP = VP(:,:,month);
                        
                        % Assume the same daily VP as the monthly mean
                        VP = repmat(VP,1,1,str2double(Tmaxfile(end-4:end-3)));
                        
                        % Get date and time info
                        time1 = ncread(Tmaxfile,'time');
                        ymd1 = char(datetime(time1/24 + datenum(1800,01,01),'ConvertFrom','datenum','Format','yyyyMMdd'))';
                        xyz.dates = ymd1;
                        xyz.time = time1;
                        
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
                        
                        % Convert from 1km to 2km resolution if required
                        if interpolate == 1
                            data = change_res(data,res2);
                        end
                        
                        % Save the derived variable
                        save_derived_nc(fname,data,xyz,Variable,template)
                    end
                end
            end
            
%         % Derived data of VP is simply the monthly mean value repeated at a
%         % daily time step for the given month, saved so that calculation of
%         % heat stress metrics in Steps 2 and 3 is more straightforward.
%         elseif strcmp(Variable,'VP')
%             % Go through each file
%             for f = 1:length(Tmaxfiles)
%                 
%                 % Get the lists of files to load for each variable
%                 Tmaxfile = [Tmaxfiles(f).folder,'/',Tmaxfiles(f).name];
%                 
%                 % Find month and year that has been loaded for daily data
%                 month = str2double(Tmaxfile(end-6:end-5));
%                 year = Tmaxfile(end-10:end-7);
%                 
%                 % 1960 is not available for VP, so skip and start at 1961
%                 if str2double(year) > 1960
%                     
%                     % Load the appropriate VP file
%                         VPfile = [VPdir,'pv_hadukgrid_uk_',res,'_mon_',year,'01-',year,'12.nc'];
%                         VP = ncread(VPfile,'pv');
%                         
%                         % Subset the correct month
%                         VP = VP(:,:,month);
%                         
%                         % Assume the same daily VP as the monthly mean
%                         data = repmat(VP,1,1,str2double(Tmaxfile(end-4:end-3)));
%                         
%                         % Get date and time info
%                         time1 = ncread(Tmaxfile,'time');
%                         ymd1 = char(datetime(time1/24 + datenum(1800,01,01),'ConvertFrom','datenum','Format','yyyyMMdd'))';
%                         xyz.dates = ymd1;
%                         xyz.time = time1;
%                     
%                         % Convert from 1km to 2km resolution if required
%                         if interpolate == 1
%                             data = change_res(data,res2);
%                         end
%                         
%                         % Save the derived variable
%                         save_derived_nc(fname,data,xyz,Variable,template)
%                 end
%             end
            
            % Otherwise, if loading a daily min. or max. heat stress variable
        else
            if strcmp(Tvar,'tasmax')
                Tfiles = Tmaxfiles;
            elseif strcmp(Tvar,'tasmin')
                Tfiles = Tminfiles;
            end
            
            % Go through each file
            for f = 1:length(Tfiles)
                
                % Get the lists of files to load for each variable
                Tfile = [Tfiles(f).folder,'/',Tfiles(f).name];
                
                % Find month and year that has been loaded for daily data
                month = str2double(Tfile(end-6:end-5));
                year = Tfile(end-10:end-7);
                
                % 1960 is not available for VP, so skip and start at 1961
                if str2double(year) > 1960
                    
                    % Load daily max and min temperatures
                    T = ncread(Tfile,Tvar);
                    time1 = ncread(Tfile,'time');
                    ymd1 = char(datetime(time1/24 + datenum(1800,01,01),'ConvertFrom','datenum','Format','yyyyMMdd'))';
                    
                    % Load the appropriate VP file
                    VPfile = [VPdir,'pv_hadukgrid_uk_',res,'_mon_',year,'01-',year,'12.nc'];
                    VP = ncread(VPfile,'pv');
                    
                    % Subset the correct month
                    VP = VP(:,:,month);
                    
                    % Assume the same daily VP as the monthly mean
                    VP = repmat(VP,1,1,length(T(1,1,:)));
                    
                    % Get date and time info
                    xyz.dates = ymd1;
                    xyz.time = time1;
                    
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
                    
                    % Save the derived variable
                    save_derived_nc(fname,data,xyz,Variable,template)
                end
            end
            
            
        end
    end
end
