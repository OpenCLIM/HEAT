function [data,xyz] = load_data(DataType,Simulation,Variable,Years,AnnSummer,TempRes)


% Domain_cell = inputs1.Domain;
% Domain_string = sprintf('%s',Domain_cell{:});

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
%         xyz.proj_x = load_UKCP_data(model,'projection_x_coordinate');
%         xyz.proj_y = load_UKCP_data(model,'projection_y_coordinate');
    else
        if strcmp(Simulation(1:2),'CP')
            model = ['cpm85',Simulation(5:6)];
            ht = ht_CPM;
%             xyz.proj_x = load_UKCP_data(model,'grid_latitude');
%             xyz.proj_y = load_UKCP_data(model,'grid_longitude');
        else
            if strcmp(Simulation(1:2),'GC')
                model = ['gcm85',Simulation(5:6)];
                ht = ht_GCM;
%                 xyz.proj_x = load_UKCP_data(model,'projection_x_coordinate');
%                 xyz.proj_y = load_UKCP_data(model,'projection_y_coordinate');
            else
                if strcmp(Simulation(1:2),'CM')
                    model = ['gcm85',Simulation(7:8)];
                    ht = ht_GCM;
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
                [data1] = load_UKCP_data(model,var);
                var = 'psl';
                [data2] = load_UKCP_data(model,var);
                data3 = p_surf(data2,data1,ht);
                
                var = 'huss';
                [data4] = load_UKCP_data(model,var);
                [data] = VapourPressure(data4,data3);
                
                
                if strcmp(Variable,'sWBGT')
                    var = 'tasmax';
                    [data5] = load_UKCP_data(model,var);
                    data = SWBGTVP(data5,data);
                    
                else
                    if strcmp(Variable,'HD')
                        var = 'tasmax';
                        [data5] = load_UKCP_data(model,var);
                        data = HumidexVP(data5,data);
                        
                    else
                        var = 'tasmax';
                        [data5] = load_UKCP_data(model,var);
                        data = AppTempVP(data5,data);
                        
                    end
                    
                end
                disp('-----')
                
            end
        end
    end
    
%     figure
%     imagesc(mean(data,3))
    
    % Load dates
    xyz.dates = load_UKCP_data(model,'yyyymmdd');
%     xyz.long = load_UKCP_data(model,'longitude');
%     xyz.lat = load_UKCP_data(model,'latitude');
    
    
    
    
end


% % Load obs data if required
% if strcmp(DataType,'obs')
%     
%     % Load each required simulation
%     for s = 1:length(inputs1.Resolution)
%         
%         % Convert cells into useable string
%         Resolution_cell = inputs1.Resolution(s);
%         Resolution_string = sprintf('%s',Resolution_cell{:});
%         
%         if ismember(inputs1.Resolution(s),'025deg')
%             disp('Error: 0.25° resolution is unavailable as an option for HadUK-Grid')
%         else
%             disp(['-> ',Resolution_string])
%         end
%     end
% end
% 
% 
% % Load obs data if required
% if ismember(inputs1.DataType(d),'reanal')
%     
%     % Load each required simulation
%     for s = 1:length(inputs1.Resolution)
%         
%         % Convert cells into useable string
%         Resolution_cell = inputs1.Resolution(s);
%         Resolution_string = sprintf('%s',Resolution_cell{:});
%         
%         if ismember(inputs1.Resolution(s),'1km')
%             disp('Error: 1km resolution is unavailable as an option for ERA5')
%         else
%             disp(['-> ',Resolution_string])
%         end
%     end
% end



