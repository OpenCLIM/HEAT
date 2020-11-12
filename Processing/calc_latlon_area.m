function [areas] = calc_latlon_area(lats,lons,method)
% [areas] = calc_latlon_area(lats,lons,method)
%
% Calculate the absolute area of each grid cell using lat-long information. 
% Two methods can be used:
%   'areaquad' to use the in built Matlab function (only suitable for
%       regular lat-long grids)
%   'geodarea' to use the downloaded geographiclib toolbox which can
%       calculate for irregular grids, but loses a line of grid cells along
%       each margin of the domain (see reference below)
% 
% Charles Karney (2019). geographiclib (https://www.mathworks.com/matlabcentral/fileexchange/50605-geographiclib), 
% MATLAB Central File Exchange. Retrieved November 26, 2019.

%% Set defaults and create output variable
if ~exist('method','var')
    disp('No method defined: deciding based on longitude grid')
    if lons(1,1) == lons(1,2) % If grid is regular
        method = 'areaquad';
    else
        method = 'geodarea';
    end
end

% Create output file
areas = nan(size(lons));


%% If using areaquad method
if strcmp(method,'areaquad')
    
    % Check if method is appropriate
    if lons(1,1) ~= lons(1,2)
        disp('Warning: Grid appears to be irregular - consider using geodarea method')
    end
    
    % Set path/defauts
    earthellipsoid = referenceSphere('earth','km');
    
    % Find grid cell edges
    for i = 1:length(lons(:,1))
        for j = 1:length(lons(1,:))
            if j < length(lons(1,:))
                lat1 = lats(i,j)+(lats(i,j)-lats(i,j+1))/2;
                lat2 = lats(i,j)-(lats(i,j)-lats(i,j+1))/2;
            else
                lat1 = lats(i,j)+(lats(i,j-1)-lats(i,j))/2;
                lat2 = lats(i,j)-(lats(i,j-1)-lats(i,j))/2;
            end
            
            if i < length(lons(:,1))
                lon1 = lons(i,j)+(lons(i,j)-lons(i+1,j))/2;
                lon2 = lons(i,j)-(lons(i,j)-lons(i+1,j))/2;
            else
                lon1 = lons(i,j)+(lons(i-1,j)-lons(i,j))/2;
                lon2 = lons(i,j)-(lons(i-1,j)-lons(i,j))/2;
            end
            
            % Calculate area
            areas(i,j) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
        end
    end
    
    
%% If using geodarea method    
else
    if strcmp(method,'geodarea')
        
        % Check if method is appropriate
        if lons(1,1) == lons(1,2)
            disp('Warning: Grid appears to be regular - consider using areaquad method')
        end
        
        % Set path/defauts
        addpath('geographiclib_toolbox_1/')
        
        % Create blank fields for corners/vertices of each grid cell
        vertices_lat = nan([size(lons),4]);
        vertices_lon = nan([size(lons),4]);
        
        % Find the corner of each gridcell
        for i = 2:length(lons(:,1))-1
            for j = 2:length(lons(1,:))-1
                
                
                vertices_lat(i,j,1) = (lats(i,j) + lats(i+1,j) + lats(i,j+1) + lats(i+1,j+1))/4;
                vertices_lat(i,j,2) = (lats(i,j) + lats(i+1,j) + lats(i,j-1) + lats(i+1,j-1))/4;
                vertices_lat(i,j,3) = (lats(i,j) + lats(i-1,j) + lats(i,j-1) + lats(i-1,j-1))/4;
                vertices_lat(i,j,4) = (lats(i,j) + lats(i-1,j) + lats(i,j+1) + lats(i-1,j+1))/4;

                vertices_lon(i,j,1) = (lons(i,j) + lons(i+1,j) + lons(i,j+1) + lons(i+1,j+1))/4;
                vertices_lon(i,j,2) = (lons(i,j) + lons(i+1,j) + lons(i,j-1) + lons(i+1,j-1))/4;
                vertices_lon(i,j,3) = (lons(i,j) + lons(i-1,j) + lons(i,j-1) + lons(i-1,j-1))/4;
                vertices_lon(i,j,4) = (lons(i,j) + lons(i-1,j) + lons(i,j+1) + lons(i-1,j+1))/4;
                
                % Calculate the area
                areas(i,j) = geodarea([vertices_lat(i,j,4) vertices_lat(i,j,3) vertices_lat(i,j,2) vertices_lat(i,j,1)], ... 
                    [vertices_lon(i,j,4) vertices_lon(i,j,3) vertices_lon(i,j,2) vertices_lon(i,j,1)]) / 1e6; % Convert to km2
            end
        end
    end
end
        
