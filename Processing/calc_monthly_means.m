function [var_months] = calc_monthly_means(var)
% [var_months] = calc_monthly_means(var)
% 
% Find monthly means of a given UKCP18 model variable to compare with
% HadUK-Grid VP or Tmean data.
% 
% Output:
%   var_months = the variable, with monthly average at daily resolution
% 
% Inputs:
%   var = the climate variable to process - this should be for the JJAS
%       summer


%% Set default year length, etc.

% Set how long the summer is for each model calendar
if rem(length(var(1,1,:)),120) == 0
    sumlen = 120;
else
    sumlen = 122;
end

% Calculate how many months there are
n_months = length(var(1,1,:))/sumlen * 4;


%% Create empty variable
var_months = nan(size(var));


%% Calculate the monthly mean
startday = 1;

% Go through each month
for i = 1:n_months
    
    % Find out how long the month is:
    % For HadGEM3 simulations with a 360 day calendar
    if sumlen == 120
        monlen = 30;
        
    else % For other calendars
        if rem(i,4) == 1 || rem(i,4) == 0 % i.e. June or September
            monlen = 30;
        else
            if rem(i,4) == 2 || rem(i,4) == 3 % i.e. July or August
                monlen = 31;
            end
        end
    end
    
    % Find start and end of month
    month_start = startday;
    month_end = startday + monlen - 1;
    
    % Set startday for next month
    startday = startday + monlen;
    
    % Calculate monthly mean
    var_month = mean(var(:,:,month_start:month_end),3);
    var_months(:,:,month_start:month_end) = repmat(var_month,1,1,monlen);


 end   

