function [data,xyz,template] = load_data(DataType,Simulation,Variable)



disp('Loading data')

%% Set directory paths
init_HEAT

% Load elevation for calculating surface pressure
generate_elevation
% Load latitude and longitude data
generate_UK_latlon
        
% Load model data if required
if strcmp(DataType,'model')
    
    % Load  required simulation
    disp(['-> Loading ',Simulation])
    
    % Change simulation name to UKCP18 file format
    % and select elevation data at correct resolution
    if strcmp(Simulation(1:2),'RC')
        model = ['rcm85',Simulation(5:6)];
        ht = ht_RCM;
        template = [UKCP18dir,'12km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-rcm_uk_12km_',model(6:7),...
                '_day_19801201-19901130.nc'];
%         xyz.proj_x = load_UKCP_data(model,'projection_x_coordinate');
%         xyz.proj_y = load_UKCP_data(model,'projection_y_coordinate');
    else
        if strcmp(Simulation(1:2),'CP')
            model = ['cpm85',Simulation(5:6)];
            ht = ht_CPM;
            template = [UKCP18dir,'2km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-cpm_uk_2.2km_',model(6:7),...
                '_day_19801201-19811130.nc'];
%             xyz.proj_x = load_UKCP_data(model,'grid_latitude');
%             xyz.proj_y = load_UKCP_data(model,'grid_longitude');
        else
            if strcmp(Simulation(1:2),'GC')
                model = ['gcm85',Simulation(5:6)];
                ht = ht_GCM;
                template = [UKCP18dir,'60km/tasmax/run',Simulation(5:6),'/tasmax_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),...
                '_day_19791201-19891130.nc'];
%                 xyz.proj_x = load_UKCP_data(model,'projection_x_coordinate');
%                 xyz.proj_y = load_UKCP_data(model,'projection_y_coordinate');
            else
                if strcmp(Simulation(1:2),'CM')
                    model = ['gcm85',Simulation(7:8)];
                    ht = ht_GCM;
                    template = [UKCP18dir,'60km/tasmax/run',Simulation(7:8),'/tasmax_rcp',model(4:5),'_land-gcm_uk_60km_',model(6:7),...
                '_day_19791201-19891130.nc'];
%                     xyz.proj_x = load_UKCP_data(model,'projection_x_coordinate');
%                     xyz.proj_y = load_UKCP_data(model,'projection_y_coordinate');
                end
            end
        end
    end
    
    % Load required variable
    disp(['---> ',Variable])
    
    if strcmp(Variable,'Tmean')
        var = 'tas';
        [data] = load_UKCP_data(model,var);
        disp('-----')
        
    else
        if strcmp(Variable,'Tmax')
            var = 'tasmax';
            [data] = load_UKCP_data(model,var);
            disp('-----')
            
        else
            if strcmp(Variable,'Tmin')
                var = 'tasmin';
                [data] = load_UKCP_data(model,var);
                disp('-----')
                
            else % All other variables will require Vapour Pressure
                
                var = 'tas';
                [data1,ymd1] = load_UKCP_data(model,var);
                var = 'psl';
                [data2,ymd2] = load_UKCP_data(model,var);
                
                [data1,ymd1,data2,ymd2] = check_consistent_timestep(data1,ymd1,data2,ymd2);
                data3 = p_surf(data2,data1,ht);
                
                var = 'huss';
                [data4,ymd4] = load_UKCP_data(model,var);
                [data3,ymd3,data4,ymd4] = check_consistent_timestep(data3,ymd1,data4,ymd4);
                [data] = VapourPressure(data4,data3);
                
                
                if strcmp(Variable,'sWBGT')
                    var = 'tasmax';
                    [data5,ymd5] = load_UKCP_data(model,var);
                    [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                    data = SWBGTVP(data5,data);
                    
                else
                    if strcmp(Variable,'HD')
                        var = 'tasmax';
                        [data5,ymd5] = load_UKCP_data(model,var);
                        [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                        data = HumidexVP(data5,data);
                        
                    else
                        var = 'tasmax';
                        [data5,ymd5] = load_UKCP_data(model,var);
                        [data,ymd,data5,ymd5] = check_consistent_timestep(data,ymd3,data5,ymd5);
                        data = AppTempVP(data5,data);
                        
                    end
                    
                end
                disp('-----')
                
            end
        end
    end
     
    % Load dates
    xyz.dates = load_UKCP_data(model,'yyyymmdd');
    xyz.time = load_UKCP_data(model,'time');

    % If variables are of inconsistent lengths, correct this
    [d1,t1,d2,t2] = check_consistent_timestep_2d(xyz.time',xyz.dates,ymd,ymd);
    xyz.dates = t1;
    xyz.time = d1;
    
    
end



