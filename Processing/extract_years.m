function [varyrs,selectedtimes] = extract_years(var,dates,startyr,endyr)
% [varyrs,selectedtimes] = extract_years(var,dates,startyr,endyr)
% 
% Take only selected years from a time series of daily data. Note: years
% may start on 1st December of the year before so they include a full
% meteorological winter at the start (i.e. DJF). This is the case for
% UKCP18 data. This function should be able to allow for this.
% 
% Output:
%   varyrs - only selected years from the input variable
%   selectedtimes - the time steps of the kept period

% Input:
%   var - the full time series of a given variable
%   dates - the time step dates for the given variable
%   startyr - the first year to be kept
%   endyr - the final year to be kept
% 

%% Find if data includes the December of the year before
if strcmp(dates(5:6,1)','12')
    yrstart = 'Dec';
else
    yrstart = 'Jan';
end


%% Find the dating index for the desired start and end
% Daily data has yyyymmdd indexing
if length(num2str(str2double(dates(:,1)))) == 8
    if strcmp(yrstart,'Dec') % UKCP18 data format starts in December
        truestart = str2double([num2str(startyr-1),'1201']);
        trueend = str2double([num2str(endyr),'1130']);
    else % WFDEI data format starts in January
        truestart = str2double([num2str(startyr),'0101']);
        trueend = str2double([num2str(endyr),'1231']);
    end
else
    % Monthly data has yyyymm indexing
    if length(num2str(str2double(dates(:,1)))) == 6
        truestart = str2double([num2str(startyr-1),'12']);
        trueend = str2double([num2str(endyr),'11']);
    end
end


%% Find when the input variable is at the start and end of desired period
for i = 1:length(dates(1,:))
    if str2double(dates(1:8,i)) == truestart
        startid = i;
    end
    
    if str2double(dates(1:8,i)) == trueend
        endid = i;
    end
end

% If not a full annual cycle is available in data, truestart and trueend
% might not yet be found. In this case, find the first date for the given
% start year and the last date for the given end year:
if ~exist('startid','var')
    startid = nan;
    i = 1;
    while isnan(startid)
        if str2double(dates(1:4,i)) == startyr
            startid = i;
        end
        i = i+1;
    end
end

if ~exist('endid','var')
    endid = nan;
    i = 1;
    while isnan(endid)
        if str2double(dates(1:4,i)) == endyr+1
            endid = i-1;
        else
            if i == length(dates(1,:))
                endid = i;
            end
        end
        i = i+1;
    end
end



%% Extract the corresponding days for the desired period
% If a climate variable has been input, subset it
if ndims(var)>2
    varyrs = var(:,:,startid:endid);
else % If not, create a blank variable
    varyrs = [];
end

% Also extract the time step info (useful for later analysis, e.g. summer removal)
selectedtimes = dates(:,startid:endid);
