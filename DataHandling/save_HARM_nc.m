function [] = save_HARM_nc(fname,data,xyz,Variable,template)
% [] = save_HARM_nc(fname,data,xyz,Variable,template)
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
% fname_long = [Deriveddir,fname,''];
fname_long = fname;

% Set some basic meta data for key derived variables
if strcmp(Variable,'Tmean')
    units = '°C';
    standard_name = 'Tmean';
    long_name = 'Daily mean temperature';
    description = 'Daily mean temperature derived from daily max. and min.';
    label_units = '°C';
    plot_label = 'Temperature at 1.5m (°C)';
end


%% Start to save
disp(['Saving workflow netCDF: ',fname])


% Load lat, long and time info for saving
x = xyz.projection_x_coordinate;
y = xyz.projection_y_coordinate;
z = xyz.dates;

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
% units = ncreadatt(template,'time','units');
% cal = ncreadatt(template,'time','calendar');

nccreate(fname_long,'time','Dimensions',{'time',length(data(1,1,:))},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'time',xyz.times);
ncwriteatt(fname_long,'time','standard_name','time');
ncwriteatt(fname_long,'time','units','hours since 1970-01-01 00:00:00');
ncwriteatt(fname_long,'time','calendar','365_day');
ncwriteatt(fname_long,'time','axis','T');

nccreate(fname_long,'yyyymmdd','Dimensions',{'time',length(data(1,1,:)),'string64',64},'Datatype','char')
ncwrite(fname_long,'yyyymmdd',xyz.dates');
ncwriteatt(fname_long,'yyyymmdd','long_name','yyyymmdd');
ncwriteatt(fname_long,'yyyymmdd','units','1');

% Write some general attributes
ncwriteatt(fname_long,'/','collection','HEAT derived variable')
ncwriteatt(fname_long,'/','creation_date',datestr(now))
%     ncwriteatt(fname_long,'/','domain',fname(11:12))
ncwriteatt(fname_long,'/','title','Variable extracted from UKCP18')
ncwriteatt(fname_long,'/','version','HEAT v1.1')




disp(['Saving of ',fname,' complete, ready for use in HARM'])
disp('-----')
