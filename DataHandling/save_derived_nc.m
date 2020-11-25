function [] = save_derived_nc(fname,data,xyz,Variable,template)

init_HEAT
% x = length(xyz.proj_x);
% y = length(xyz.proj_x);
% z = length(xyz.dates(1,:));


% Set some metadata
if strcmp(Variable,'VP')
    units = 'mb';
    standard_name = 'vapour_pressure';
    long_name = 'Vapour pressure';
    description = 'Vapour Pressure';
    label_units = 'mb';
    plot_label = 'Vapour pressure at 1.5m (mb)';
end

if strcmp(Variable,'AT')
    units = '°C';
    standard_name = 'apparent_temperature';
    long_name = 'Apparent temperature';
    description = 'Apparent temperature';
    label_units = '°C';
    plot_label = 'Apparent temperature at 1.5m (°C)';
end

if strcmp(Variable,'HD')
    units = 'N/A';
    standard_name = 'humidex';
    long_name = 'Humidex';
    description = 'Humidex';
    label_units = 'N/A';
    plot_label = 'Humidex at 1.5m (dimensionless)';
end

if strcmp(Variable,'sWBGT')
    units = '°C';
    standard_name = 'sWBGT';
    long_name = 'Simplified Wet Buld Globe Temperature';
    description = 'Simplified Wet Buld Globe Temperature';
    label_units = '°C';
    plot_label = 'sWBGT at 1.5m (mb)';
end


disp(['Saving derived variable: ',fname])

template_data = ncinfo(template);

x = ncread(template,'projection_x_coordinate');
y = ncread(template,'projection_y_coordinate');
z = ncread(template,'time');


nccreate([Deriveddir,fname],Variable,'Dimensions',{'projection_x_coordinate',length(x),'projection_y_coordinate',length(y),'time',length(data(1,1,:))},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite([Deriveddir,fname],Variable,data);
ncwriteatt([Deriveddir,fname],Variable,'standard_name',standard_name);
ncwriteatt([Deriveddir,fname],Variable,'long_name',long_name);
ncwriteatt([Deriveddir,fname],Variable,'units',units);
ncwriteatt([Deriveddir,fname],Variable,'description',description);
ncwriteatt([Deriveddir,fname],Variable,'label_units',label_units);
ncwriteatt([Deriveddir,fname],Variable,'plot_label',plot_label);
ncwriteatt([Deriveddir,fname],Variable,'grid_mapping','transverse_mercator');
ncwriteatt([Deriveddir,fname],Variable,'coordinates','ensemble_member_id latitude longitude month_number year yyyymmdd');
ncwriteatt([Deriveddir,fname],Variable,'missing_value',-9999);

nccreate([Deriveddir,fname],'projection_x_coordinate','Dimensions',{'projection_x_coordinate',length(x)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
ncwrite([Deriveddir,fname],'projection_x_coordinate',x);
ncwriteatt([Deriveddir,fname],'projection_x_coordinate','standard_name','easting');
ncwriteatt([Deriveddir,fname],'projection_x_coordinate','long_name','easting');
ncwriteatt([Deriveddir,fname],'projection_x_coordinate','units','m');
ncwriteatt([Deriveddir,fname],'projection_x_coordinate','axis','X');

nccreate([Deriveddir,fname],'projection_y_coordinate','Dimensions',{'projection_y_coordinate',length(y)},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
ncwrite([Deriveddir,fname],'projection_y_coordinate',y);
ncwriteatt([Deriveddir,fname],'projection_y_coordinate','standard_name','northing');
ncwriteatt([Deriveddir,fname],'projection_y_coordinate','long_name','northing');
ncwriteatt([Deriveddir,fname],'projection_y_coordinate','units','m');
ncwriteatt([Deriveddir,fname],'projection_y_coordinate','axis','Y');

nccreate([Deriveddir,fname],'time','Dimensions',{'time',length(data(1,1,:))},'Datatype','single','Format','netcdf4_classic','DeflateLevel',2)
ncwrite([Deriveddir,fname],'time',xyz.time);
ncwriteatt([Deriveddir,fname],'time','standard_name','time');
ncwriteatt([Deriveddir,fname],'time','units','hours since 1970-01-01 00:00:00');
ncwriteatt([Deriveddir,fname],'time','calendar','360_day');
ncwriteatt([Deriveddir,fname],'time','axis','T');


nccreate([Deriveddir,fname],'yyyymmdd','Dimensions',{'time',length(data(1,1,:)),'string64',64},'Datatype','char')
ncwrite([Deriveddir,fname],'yyyymmdd',xyz.dates');
ncwriteatt([Deriveddir,fname],'yyyymmdd','long_name','yyyymmdd');
ncwriteatt([Deriveddir,fname],'yyyymmdd','units','1');

ncwriteatt([Deriveddir,fname],'/','collection','HEAT derived variable')
ncwriteatt([Deriveddir,fname],'/','creation_date',datestr(now))
ncwriteatt([Deriveddir,fname],'/','domain',fname(11:12))
ncwriteatt([Deriveddir,fname],'/','title','Vairable derived from UKCP18')
ncwriteatt([Deriveddir,fname],'/','version',fname(6:9))

% copyfile(template,[Deriveddir,'/',fname])
%
%
% ncid = netcdf.open([Deriveddir,'/',fname],'NC_WRITE');
%
% netcdf.reDef(ncid)
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,0)
% netcdf.renameVar(ncid,0,standard_name)
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,0)
%
% % netcdf.putVar(ncid,varid,data)
%
% netcdf.close(ncid)

disp(['Loading, processing and saving of ',Variable,' complete'])
disp('-----')
