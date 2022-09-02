function [data,dates,times] = subset_temporal(data,dates,times,Years,AnnSummer)


%% Subset the correct temporal period
% Subset the correct years
tsstart = find_date(num2str(Years(1)),dates(1:8,:));
tsend = find_date(num2str(Years(2)),dates(1:8,:));
data = data(:,:,tsstart:tsend);
dates = dates(:,tsstart:tsend);
times = times(tsstart:tsend);

% Subset summers if necessary
if strcmp(AnnSummer,'Summer')
    summer = 'MO';
%     disp('Subsetting summer (1st June-15th Sept.) days')
%     disp('-----')
    data = extract_summers(data,dates,summer);
    dates = extract_summers(dates,dates,summer);
    times = extract_summers(times,dates,summer);
else
    if strcmp(AnnSummer,'JJA')
        summer = 'JJA';
%         disp('Subsetting JJA days')
%         disp('-----')
        data = extract_summers(data,dates,summer);
        dates = extract_summers(dates,dates,summer);
        times = extract_summers(times,dates,summer);

    end
end


