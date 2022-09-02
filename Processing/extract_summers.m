function [varsum] = extract_summers(var,dates,summertype)
% [varsum] = extract_summers(var,dates,summertype)
% 
% Take only the summer days (either JJA or Met. Office definition) from a
% time series. If only a selected time period (e.g. 1981-2010) is wanted,
% this should be done with extract_years.m before extracting the summer
% days.
% 
% Output:
%   varsum - only summer days from the input variable

% Input:
%   var - the full time series of a given variable
%   dates - the 'yyyymmdd' data for the input variable (if data has already
%       been trimmed to a certain time period using extract_years.m, this
%       should be the selectedtimes output variable)
%   summertype - the desired summer period ('JJA' = June July August, 'MO'
%       = Met. Office summer definition (1st June-15th September), 'JJAS' =
%       June July August September)
% 

%% Find the summer days of interest

% Generate empty array to store the indexes for summer days
summerid = [];

% Both summers include all of June, July and August
for i = 1:length(dates(1,:))
    if strcmp(num2str(str2double(dates(5:6,i))),'6')
        summerid = cat(1,summerid,i);
    else
        if strcmp(num2str(str2double(dates(5:6,i))),'7')
            summerid = cat(1,summerid,i);
        else
            if strcmp(num2str(str2double(dates(5:6,i))),'8')
                summerid = cat(1,summerid,i);
            end
        end
    end
    
    % Find September days if using the Met. Office summer definition
    if strcmp(summertype,'MO')
        if strcmp(num2str(str2double(dates(5:8,i))),'901')
            summerid = cat(1,summerid,[i:i+14]');
        end
    end
    
    % Find September days if using the JJAS summer definition
    if strcmp(summertype,'JJAS')
        if strcmp(num2str(str2double(dates(5:6,i))),'9')
            summerid = cat(1,summerid,i);
        end
    end
end


%% Extract the corresponding days for the desired period
if ndims(var) == 3
    varsum = var(:,:,summerid);
elseif ndims(var) == 2
    varsum = var(:,summerid);
elseif ndims(var) == 1
    varsum = var(summerid);
end



