function [data,xyz,template] = load_data(DataType,Simulation,Variable,CPM_period)
% [data,xyz,template] = load_data(DataType,Simulation,Variable,CPM_period)
% 
% Load the raw data for HEAT calculation of heat stress metrics and saving
% derived data.
% 
% Inputs:
%  - DataType = 'model', 'obs' or 'reanal' as defined in the input file
%  - Simulation = simulation name (e.g. RCM-01) if using model simulations,
%   or the resolution if using observations or reanalysis data.
%  - Variable = the climate variable to load
%  - CPM_period = only required if loading CPM data, in which case the
%   20-year period to be loaded can be specified (default = 1981-2000)
% 
% Outputs:
%  - data = the loaded dataset
%  - xyz = some lat-long and time data about the loaded data
%  - template = a charcter string directing to a variable with a similar
%   netCDF format, which is used when saving derived variables

disp('Loading data')

%% Set directory paths
init_HEAT

% Load elevation and lat-long data
load_xyz
        

%% Load UKCP18 data if required
if strcmp(DataType,'UKCP18')
    
    % Load  required simulation
    disp(['-> Loading ',Simulation])
        
    % Change simulation name to UKCP18 file format
    % and select elevation data at correct resolution
    if strcmp(Simulation(1:2),'RC')
        
        % Find raw data location
        model = ['rcm85',Simulation(5:6)];
        ht = ht_12km;
        Rawdir = [UKCP18dir,'12km/tasmax/run',Simulation(5:6),'/'];
        
        template = [UKCP18dir,'12km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),...
            '_day_19801201-19901130.nc'];
        
    elseif strcmp(Simulation(1:2),'CP')
        model = ['cpm85',Simulation(5:6)];
        ht = ht_2km;
        Rawdir = [UKCP18dir,'2km/tasmax/run',Simulation(5:6),'/'];
        
        template = [UKCP18dir,'2km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-cpm_uk_2.2km_',model(6:7),...
            '_day_19801201-19811130.nc'];
        
    elseif strcmp(Simulation(1:2),'GC')
        model = ['gcm85',Simulation(5:6)];
        ht = ht_60km;
        Rawdir = [UKCP18dir,'60km/tasmax/run',Simulation(5:6),'/'];
        
        template = [UKCP18dir,'60km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),...
            '_day_19791201-19891130.nc'];
        
    elseif strcmp(Simulation(1:2),'CM')
        model = ['gcm85',Simulation(7:8)];
        ht = ht_60km;
        Rawdir = [UKCP18dir,'60km/tasmax/run',Simulation(5:6),'/'];
        
        template = [UKCP18dir,'60km/tasmax/run',Simulation(7:8),'/tasmax_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),...
            '_day_19791201-19891130.nc'];
        
    end

    end
    
    % Find how many files are to be loaded/produced
    files = dir([Rawdir '*.nc']);
    
    for f = 1:length(files)
        
        file = [files(f).folder,'/',files(f).name];
    
    % Load required variable
    disp(['---> ',Variable])
    
    if strcmp(Variable,'Tmean')
        var = 'tas';
        [data] = load_UKCP_data(model,var,CPM_period);
        disp('-----')
        
    else
        if strcmp(Variable,'Tmax')
            var = 'tasmax';
            [data] = load_UKCP_data(model,var,CPM_period);
            disp('-----')
            
        else
            if strcmp(Variable,'Tmin')
                var = 'tasmin';
                [data] = load_UKCP_data(model,var,CPM_period);
                disp('-----')
                
            else % All other variables will require Vapour Pressure
                
                var = 'tas';
                [data1,ymd1] = load_UKCP_data(model,var,CPM_period);
                var = 'psl';
                [data2,ymd2] = load_UKCP_data(model,var,CPM_period);
                
                [data1,ymd1,data2,ymd2] = check_consistent_timestep(data1,ymd1,data2,ymd2);
                data3 = p_surf(data2,data1,ht);
                
                var = 'huss';
                [data4,ymd4] = load_UKCP_data(model,var,CPM_period);
                [data3,ymd3,data4,ymd4] = check_consistent_timestep(data3,ymd1,data4,ymd4);
                [data] = VapourPressure(data4,data3);
                
                
                if strcmp(Variable,'sWBGT')
                    var = 'tasmax';
                    [data5,ymd5] = load_UKCP_data(model,var,CPM_period);
                    [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                    data = SWBGTVP(data5,data);
                    
                else
                    if strcmp(Variable,'HD')
                        var = 'tasmax';
                        [data5,ymd5] = load_UKCP_data(model,var,CPM_period);
                        [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                        data = HumidexVP(data5,data);
                        
                    else
                        var = 'tasmax';
                        [data5,ymd5] = load_UKCP_data(model,var,CPM_period);
                        [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                        data = AppTempVP(data5,data);
                        
                    end
                    
                end
                disp('-----')
                
            end
        end
    end
     
    % Load dates
    xyz.dates = load_UKCP_data(model,'yyyymmdd',CPM_period);
    xyz.time = load_UKCP_data(model,'time',CPM_period);

    % If multiple variables have been loaded, if they were of inconsistent 
    % lengths and have been corrected, then their date strings will be
    % inconsistent ? correct this:
    if exist('ymd','var')
        [d1,t1,d2,t2] = check_consistent_timestep_2d(xyz.time',xyz.dates,ymd,ymd);
        xyz.dates = t1;
        xyz.time = d1;
    end    
end


%% If loading observational data
if strcmp(DataType,'HadUKGrid')
    
    % Load  required simulation
    disp(['-> Loading HadUK-Grid observations at ',Simulation,' resolution'])
    
    % Use correct elevation dataset for this resolution
    if strcmp(Simulation,'12km')
        ht = ht_12km;
    else
        if strcmp(Simulation,'60km')
            ht = ht_60km;
        else
            ht = ht_1km;
        end
    end
    
    % Load required variable
    disp(['---> ',Variable])
    
    if strcmp(Variable,'Tmean')
        var = 'tas';
        [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
        data = cat(3,data80s,data90s,data00s,data10s);
        disp('-----')
        
    else
        if strcmp(Variable,'Tmax')
            var = 'tasmax';
            [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
            data = cat(3,data80s,data90s,data00s,data10s);
            disp('-----')
            
        else
            if strcmp(Variable,'Tmin')
                var = 'tasmin';
                [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                data = cat(3,data80s,data90s,data00s,data10s);
                disp('-----')
                
            else % All other variables will require Vapour Pressure
                
                var = 'tas';
                [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                data1 = cat(3,data80s,data90s,data00s,data10s);
                var = 'psl';
                [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                data2 = cat(3,data80s,data90s,data00s,data10s);
                
                data3 = p_surf(data2,data1,ht);
                
                var = 'huss';
                [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                data4 = cat(3,data80s,data90s,data00s,data10s);
                [data] = VapourPressure(data4,data3);
                
                
                if strcmp(Variable,'sWBGT')
                    var = 'tasmax';
                    [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                    data5 = cat(3,data80s,data90s,data00s,data10s);
                    data = SWBGTVP(data5,data);
                    
                else
                    if strcmp(Variable,'HD')
                        var = 'tasmax';
                        [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                        data5 = cat(3,data80s,data90s,data00s,data10s);
                        data = HumidexVP(data5,data);
                        
                    else
                        var = 'tasmax';
                        [data80s,data90s,data00s,data10s] = load_HadUK_data(var,Simulation);
                        data5 = cat(3,data80s,data90s,data00s,data10s);
                        data = AppTempVP(data5,data);
                        
                    end
                    
                end
                disp('-----')
                
            end
        end
    end
     
    % Load time info from netCDF and convert to yyyymmdd 
    xyz.time = load_HadUK_data('time',Simulation);
    xyz.dates = char(datetime(xyz.time/24 + 657438,'ConvertFrom','datenum','Format','yyyyMMdd')); % 657438 is the offset to the HadUK-Grid start date, obtained by: datenum([1800 01 01 00 00 00])
    xyz.dates = xyz.dates';
    
    % Provide directory path to appropriate template file
%     if strcmp(Simulation(1:2),'12')
        template = [HadUKdir,'v1.0.2.1/tasmax/',Simulation,'/tasmax_hadukgrid_uk_',Simulation,'_day_19810501-19810531.nc'];
%     else
%         if strcmp(Simulation(1:2),'60')
%         template = [HadUKdir,'v1.0.2.1/tasmax/',Simulation,'/tasmax_hadukgrid_uk_',Simulation,'_day_19810501-19810531.nc'];
%     end
end



