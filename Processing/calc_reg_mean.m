function [regmean] = calc_reg_mean(data,reg)
% Calculate the regional mean of dataset 'data' for region 'reg'



% Regions for reference:
regs = {'Scotland','North East','North West','Yorkshire and the Humber','East Midlands','West Midlands','East of England','Greater London','South East','South West','Wales','Northern Ireland','Isle of Mann'};

% Load the regional area data
generate_region_latlon_area

% Find what resolution of data is being used
if length(data(:,1,1)) == 17
    area = areas_GCM_frac_regions;
else
    if length(data(:,1,1)) == 82
        area = areas_RCM_frac_regions;
    else
        if length(data(:,1,1)) == 484
            area = areas_CPM_frac_regions;
        else
            if length(data(:,1,1)) == 121
                area = areas_ERA5_frac_regions;
            end
        end
    end
end

% Calculate the regional mean for the selected region
for i = 1:length(regs)
    if strcmp(regs(i),reg)
        regmean = squeeze(nansum(nansum(data .* area(:,:,i),1),2));
    end
end
