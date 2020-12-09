function [data,dates] = subset_temporal(data,dates,Years,AnnSummer)


%% Subset the correct temporal period
% Subset the correct years
[data,dates] = extract_years(data,dates',Years(1),Years(2));

% Subset summers if necessary
if strcmp(AnnSummer,'Summer')
    summer = 'MO';
    disp('Subsetting summer (1st June-15th Sept.) days')
    disp('-----')
    data = extract_summers(data,dates,summer);
    dates = extract_summers(dates,dates,summer);
else
    if strcmp(AnnSummer,'JJA')
        summer = 'JJA';
        disp('Subsetting JJA days')
        disp('-----')
        dates = extract_summers(dates,dates,summer);
    end
end


