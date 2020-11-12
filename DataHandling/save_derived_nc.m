function [] = save_derived_nc(fname,data,xyz,Variable)

% x = length(xyz.proj_x);
% y = length(xyz.proj_x);
% z = length(xyz.dates(1,:));
    
    disp(['Saving derived variable: ',fname])
    disp(['Loading, processing and saving of ',Variable,' complete'])
    disp('-----')
%     nccreate(fname,Variable,'Dimensions',{'projection_x_coordinate',x,'projection_y_coordinate',y,'yyyymmdd',z},'Datatype','double')
%     ncwrite(fname,Variable,data);
%     ncwriteatt(fname,Variable,'long_name','Recent past daily min. 2m temperature on days exceeding summer Tmax 95th percentile');
%     ncwriteatt(fname,Variable,'units','deg. C');
%     ncwriteatt(fname,Variable,'missing_value',-32767);
% 
%     
%     
%     nccreate(fname,'projection_x_coordinate','Dimensions',{'projection_x_coordinate',x},'Datatype','single')
%     ncwrite(fname,'projection_x_coordinate',xs);
%     ncwriteatt(fname,'projection_x_coordinate','standard_name','easting');
%     ncwriteatt(fname,'projection_x_coordinate','long_name','easting');
%     ncwriteatt(fname,'projection_x_coordinate','units','m');
%     ncwriteatt(fname,'projection_x_coordinate','axis','X');
%     
%     nccreate(fname,'projection_y_coordinate','Dimensions',{'projection_y_coordinate',y},'Datatype','single')
%     ncwrite(fname,'projection_y_coordinate',ys);
%     ncwriteatt(fname,'projection_y_coordinate','standard_name','northing');
%     ncwriteatt(fname,'projection_y_coordinate','long_name','northing');
%     ncwriteatt(fname,'projection_y_coordinate','units','m');
%     ncwriteatt(fname,'projection_y_coordinate','axis','Y');
    