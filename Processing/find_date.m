function [date_id] = find_date(date_str,dates)
% [date_id] = find_date(date_str)
% 
% Find the index in the time domain of a specified date.

for i = 1:length(dates(1,:))
    if strcmp(dates(:,i)',date_str)
        date_id = i;
    end
end