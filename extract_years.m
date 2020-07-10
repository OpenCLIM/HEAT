function [varyrs,selectedtimes] = extract_years(var,dates,startyr,endyr,yrstart)
% [varyrs,selectedtimes] = extract_years(var,dates,startyr,endyr,yrstart)
% 
% Take only selected years from a time series of daily data. Note: years
% are taken as starting on 1st December of the year before so they include
% a full meteorological winter at the start (i.e. DJF).
% 
% Output:
%   varyrs - only selected years from the input variable
%   selectedtimes - the time steps of the kept period

% Input:
%   var - the full time series of a given variable
%   dates - the time step dates for the given variable
%   startyr - the first year to be kept
%   endyr - the final year to be kept
%   yrstart - the start point for the year ('Dec': 1st December for UKCP18
%       simulations [default], 'Jan': 1st January for WFDEI data)
% 

%% Set default yrstart
if ~exist('yrstart','var')
    yrstart = 'Dec';
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


%% Find when the input variable is at the start and end
for i = 1:length(dates(1,:))
    if str2double(dates(1:8,i)) == truestart
        startid = i;
    end
    
    if str2double(dates(1:8,i)) == trueend
        endid = i;
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
