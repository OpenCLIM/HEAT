function [] = save_derived_nc(fname,data,xyz,Variable,template)
% [] = save_derived_nc(fname,data,xyz,Variable,template)
% 
% This function saves a 'derived variable' from HEAT, i.e. something that
% is calculated from more than one raw variable but remains on the same
% temporal and spatial domain, e.g. sWBGT. Unlike the raw UKCP18 data, the
% derived data is concatenated into a single time series. There may be some
% changes to the meta data as only the attributes that were essential to
% saving were included.
% 
% 

%% Start by loading the required data


%% Set some basic details for saving
% Load the directory paths to know where to save derived data
init_HEAT

% Find the start and end of the record
period_start = xyz.dates(1:8,1);
period_end = xyz.dates(1:8,length(xyz.dates(1,:)));

% Use the start/end info to create the file name for the derived variable
fname = [fname,'-',period_start','-',period_end','.nc'];
fname_long = [Deriveddir,fname,''];

% Set some basic meta data for key derived variables
if strcmp(Variable,'T')
    units = '°C';
    standard_name = 'Tmean';
    long_name = 'Daily mean temperature';
    description = 'Daily mean temperature derived from daily max. and min.';
    label_units = '°C';
    plot_label = 'Temperature at 1.5m (°C)';
end

if strcmp(Variable,'VP')
    units = 'mb';
    standard_name = 'vapour_pressure';
    long_name = 'Vapour pressure daily mean';
    description = 'Vapour Pressure';
    label_units = 'mb';
    plot_label = 'Vapour pressure at 1.5m (mb)';
end

if strcmp(Variable,'AT')
    units = '°C';
    standard_name = 'apparent_temperature';
    long_name = 'Apparent temperature daily mean';
    description = 'Apparent temperature';
    label_units = '°C';
    plot_label = 'Apparent temperature at 1.5m (°C)';
end

if strcmp(Variable,'HD')
    units = 'N/A';
    standard_name = 'humidex';
    long_name = 'Humidex daily mean';
    description = 'Humidex';
    label_units = 'N/A';
    plot_label = 'Humidex at 1.5m (dimensionless)';
end

if strcmp(Variable,'sWBGT')
    units = '°C';
    standard_name = 'sWBGT';
    long_name = 'Simplified Wet Buld Globe Temperature daily mean';
    description = 'Simplified Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'sWBGT at 1.5m (°C)';
end

if strcmp(Variable,'WBGT')
    units = '°C';
    standard_name = 'WBGT';
    long_name = 'Wet Buld Globe Temperature daily mean';
    description = 'Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'WBGT at 1.5m (°C)';
end


% Set some basic meta data for key derived variables
if strcmp(Variable,'VPmax')
    units = 'mb';
    standard_name = 'vapour_pressure_max';
    long_name = 'Vapour pressure daily maximum';
    description = 'Vapour Pressure';
    label_units = 'mb';
    plot_label = 'Vapour pressure at 1.5m (mb)';
end

if strcmp(Variable,'ATmax')
    units = '°C';
    standard_name = 'apparent_temperature_max';
    long_name = 'Apparent temperature daily maximum';
    description = 'Apparent temperature';
    label_units = '°C';
    plot_label = 'Apparent temperature at 1.5m (°C)';
end

if strcmp(Variable,'HDmax')
    units = 'N/A';
    standard_name = 'humidex_max';
    long_name = 'Humidex daily maximum';
    description = 'Humidex';
    label_units = 'N/A';
    plot_label = 'Humidex at 1.5m (dimensionless)';
end

if strcmp(Variable,'sWBGTmax')
    units = '°C';
    standard_name = 'sWBGT_max';
    long_name = 'Simplified Wet Buld Globe Temperature daily maximum';
    description = 'Simplified Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'sWBGT at 1.5m (°C)';
end

if strcmp(Variable,'sWBGTmax')
    units = '°C';
    standard_name = 'WBGT_max';
    long_name = 'Wet Buld Globe Temperature daily maximum';
    description = 'Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'WBGT at 1.5m (°C)';
end

% Set some basic meta data for key derived variables
if strcmp(Variable,'VPmin')
    units = 'mb';
    standard_name = 'vapour_pressure_min';
    long_name = 'Vapour pressure daily minimum';
    description = 'Vapour Pressure';
    label_units = 'mb';
    plot_label = 'Vapour pressure at 1.5m (mb)';
end

if strcmp(Variable,'ATmin')
    units = '°C';
    standard_name = 'apparent_temperature_min';
    long_name = 'Apparent temperature daily minimum';
    description = 'Apparent temperature';
    label_units = '°C';
    plot_label = 'Apparent temperature at 1.5m (°C)';
end

if strcmp(Variable,'HDmin')
    units = 'N/A';
    standard_name = 'humidex_min';
    long_name = 'Humidex daily minimum';
    description = 'Humidex';
    label_units = 'N/A';
    plot_label = 'Humidex at 1.5m (dimensionless)';
end

if strcmp(Variable,'sWBGTmin')
    units = '°C';
    standard_name = 'sWBGT_min';
    long_name = 'Simplified Wet Buld Globe Temperature daily minimum';
    description = 'Simplified Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'sWBGT at 1.5m (°C)';
end

if strcmp(Variable,'WBGTmin')
    units = '°C';
    standard_name = 'WBGT_min';
    long_name = 'Wet Buld Globe Temperature daily minimum';
    description = 'Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'WBGT at 1.5m (°C)';
end


%% Start to save
disp(['Saving derived variable: ',fname])


% template_data = ncinfo(template);

% Find if saving CPM simulation or any other kind of simulation (lat long
% is different for each)
% If data is CPM:
if length(data(:,1,1)) == 484
    
    template2 = [UKCP18dir,'2km/tasmax/run01/tasmax_rcp85_land-cpm_uk_2.2km_01_day_19801201-19811130.nc'];
    % Load lat, long and time info for saving
    y = ncread(template2,'grid_latitude');
    x = ncread(template2,'grid_longitude');
    yy = ncread(template2,'latitude');
    xx = ncread(template2,'longitude');
        
    % Create netCDF and derived variable
    nccreate(fname_long,Variable,'Dimensions',{'grid_longitude',length(x),'grid_latitude',length(y),'time',length(data(1,1,:))},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,Variable,data);
    ncwriteatt(fname_long,Variable,'standard_name',standard_name);
    ncwriteatt(fname_long,Variable,'long_name',long_name);
    ncwriteatt(fname_long,Variable,'units',units);
    ncwriteatt(fname_long,Variable,'description',description);
    ncwriteatt(fname_long,Variable,'label_units',label_units);
    ncwriteatt(fname_long,Variable,'plot_label',plot_label);
    ncwriteatt(fname_long,Variable,'grid_mapping','rotated_latitude_longitude');
    ncwriteatt(fname_long,Variable,'coordinates','ensemble_member_id latitude longitude month_number year yyyymmdd');
    ncwriteatt(fname_long,Variable,'missing_value',-9999);
    
    % Add lat and long data
    nccreate(fname_long,'latitude','Dimensions',{'grid_longitude',length(x),'grid_latitude',length(y)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'latitude',yy);
    ncwriteatt(fname_long,'latitude','standard_name','latitude');
    ncwriteatt(fname_long,'latitude','units','degrees_north');
    
    nccreate(fname_long,'longitude','Dimensions',{'grid_longitude',length(x),'grid_latitude',length(y)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'longitude',xx);
    ncwriteatt(fname_long,'longitude','standard_name','longitude');
    ncwriteatt(fname_long,'longitude','units','degrees_east');
    
    nccreate(fname_long,'grid_latitude','Dimensions',{'grid_latitude',length(y)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'grid_latitude',y);
    ncwriteatt(fname_long,'grid_latitude','standard_name','grid_latitude');
    ncwriteatt(fname_long,'grid_latitude','long_name','grid_latitude');
    ncwriteatt(fname_long,'grid_latitude','units','degrees');
    ncwriteatt(fname_long,'grid_latitude','axis','Y');
    
    nccreate(fname_long,'grid_longitude','Dimensions',{'grid_longitude',length(x)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'grid_longitude',x);
    ncwriteatt(fname_long,'grid_longitude','standard_name','grid_longitude');
    ncwriteatt(fname_long,'grid_longitude','long_name','grid_longitude');
    ncwriteatt(fname_long,'grid_longitude','units','degrees');
    ncwriteatt(fname_long,'grid_longitude','axis','X');
    
    % Load info about lat long rotations
    nccreate(fname_long,'rotated_latitude_longitude','Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwriteatt(fname_long,'rotated_latitude_longitude','grid_north_pole_latitude',37.5);
    ncwriteatt(fname_long,'rotated_latitude_longitude','grid_north_pole_longitude',177.5);
    ncwriteatt(fname_long,'rotated_latitude_longitude','north_pole_grid_longitude',0.0);
    ncwriteatt(fname_long,'rotated_latitude_longitude','earth_radius',6371229.0);
    ncwriteatt(fname_long,'rotated_latitude_longitude','grid_mapping_name','rotated_latitude_longitude');
    ncwriteatt(fname_long,'rotated_latitude_longitude','longitude_of_prime_meridian',0.0);
    
    % Add time and date data
    % Some attributes are easier to read from template than to define
    % manually as they change between simulations
    units = ncreadatt(template,'time','units');
    cal = ncreadatt(template,'time','calendar');
    
    nccreate(fname_long,'time','Dimensions',{'time',length(data(1,1,:))},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'time',xyz.time);
    ncwriteatt(fname_long,'time','standard_name','time');
    ncwriteatt(fname_long,'time','units',units);
    ncwriteatt(fname_long,'time','calendar',cal);
    ncwriteatt(fname_long,'time','axis','T');
    
    nccreate(fname_long,'yyyymmdd','Dimensions',{'time',length(data(1,1,:)),'string64',64},'Datatype','char')
    ncwrite(fname_long,'yyyymmdd',xyz.dates');
    ncwriteatt(fname_long,'yyyymmdd','long_name','yyyymmdd');
    ncwriteatt(fname_long,'yyyymmdd','units','1');
    
    % Write some general attributes
    ncwriteatt(fname_long,'/','collection','HEAT derived variable')
    ncwriteatt(fname_long,'/','creation_date',datestr(now))
%     ncwriteatt(fname_long,'/','domain',fname(11:12))
    ncwriteatt(fname_long,'/','title','Variable derived from UKCP18')
    ncwriteatt(fname_long,'/','version','HEAT v1.0')
    
    
    
  
% Otherswise follow procedure for other GCMs/RCMs
else
    
    % Load lat, long and time info for saving
    x = ncread(template,'projection_x_coordinate');
    y = ncread(template,'projection_y_coordinate');
    z = ncread(template,'time');
    
    % Create netCDF and derived variable
    nccreate(fname_long,Variable,'Dimensions',{'projection_x_coordinate',length(x),'projection_y_coordinate',length(y),'time',length(data(1,1,:))},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,Variable,data);
    ncwriteatt(fname_long,Variable,'standard_name',standard_name);
    ncwriteatt(fname_long,Variable,'long_name',long_name);
    ncwriteatt(fname_long,Variable,'units',units);
    ncwriteatt(fname_long,Variable,'description',description);
    ncwriteatt(fname_long,Variable,'label_units',label_units);
    ncwriteatt(fname_long,Variable,'plot_label',plot_label);
    ncwriteatt(fname_long,Variable,'grid_mapping','transverse_mercator');
    ncwriteatt(fname_long,Variable,'coordinates','ensemble_member_id latitude longitude month_number year yyyymmdd');
    ncwriteatt(fname_long,Variable,'missing_value',-9999);
    
    % Add lat and long data
    nccreate(fname_long,'projection_x_coordinate','Dimensions',{'projection_x_coordinate',length(x)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_x_coordinate',x);
    ncwriteatt(fname_long,'projection_x_coordinate','standard_name','easting');
    ncwriteatt(fname_long,'projection_x_coordinate','long_name','easting');
    ncwriteatt(fname_long,'projection_x_coordinate','units','m');
    ncwriteatt(fname_long,'projection_x_coordinate','axis','X');
    
    nccreate(fname_long,'projection_y_coordinate','Dimensions',{'projection_y_coordinate',length(y)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'projection_y_coordinate',y);
    ncwriteatt(fname_long,'projection_y_coordinate','standard_name','northing');
    ncwriteatt(fname_long,'projection_y_coordinate','long_name','northing');
    ncwriteatt(fname_long,'projection_y_coordinate','units','m');
    ncwriteatt(fname_long,'projection_y_coordinate','axis','Y');
    
    % Add time and date data
    % Some attributes are easier to read from template than to define
    % manually as they change between simulations
    units = ncreadatt(template,'time','units');
    cal = ncreadatt(template,'time','calendar');
    
    nccreate(fname_long,'time','Dimensions',{'time',length(data(1,1,:))},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
    ncwrite(fname_long,'time',xyz.time);
    ncwriteatt(fname_long,'time','standard_name','time');
    ncwriteatt(fname_long,'time','units',units);
    ncwriteatt(fname_long,'time','calendar',cal);
    ncwriteatt(fname_long,'time','axis','T');
    
    nccreate(fname_long,'yyyymmdd','Dimensions',{'time',length(data(1,1,:)),'string64',64},'Datatype','char')
    ncwrite(fname_long,'yyyymmdd',xyz.dates');
    ncwriteatt(fname_long,'yyyymmdd','long_name','yyyymmdd');
    ncwriteatt(fname_long,'yyyymmdd','units','1');
    
    % Write some general attributes
    ncwriteatt(fname_long,'/','collection','HEAT derived variable')
    ncwriteatt(fname_long,'/','creation_date',datestr(now))
%     ncwriteatt(fname_long,'/','domain',fname(11:12))
    ncwriteatt(fname_long,'/','title','Variable derived from UKCP18')
    ncwriteatt(fname_long,'/','version','HEAT v1.0')
    

end


disp(['Loading, processing and saving of ',fname,' complete'])
disp('-----')
