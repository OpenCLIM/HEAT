function [lon_id,lat_id] = find_location(lon_x,lat_x,lons,lats)
% [lon_id,lat_id] = find_location(lon_x,lat_x,lons,lats)
% 
% Find the lat-long index for a given location (lon_x,lat_x).
% 
% For example, London ~ lon_x = 0; lat_x = 51.5;

lons_dif = lons-lon_x;
lats_dif = lats-lat_x;
comb_dif = sqrt(lons_dif.^2+lats_dif.^2);
min_val = min(min(comb_dif));

for i = 1:length(lons_dif(:,1))
    for j = 1:length(lons_dif(1,:))
        if comb_dif(i,j) == min_val
            lon_id = i;
            lat_id = j;
        end
    end
end
