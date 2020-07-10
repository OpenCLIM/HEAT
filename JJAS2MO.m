function [MOdays] = JJAS2MO(JJASdays)
% [MOdays] = JJAS2MO(JJASdays)
% 
% Convert a time series of JJAS days to MO summer days.
% 

if rem(length(JJASdays(1,1,:)),120) == 0
    % Model is using HadGEM3 360-day calendar
    
    % Find how many years are in time series and make empty array
    yrs = length(JJASdays(1,1,:))/120;
    MOdays = zeros(length(JJASdays(:,1,1)),length(JJASdays(1,:,1)),yrs*105);
    
    for i = 1:yrs
        
        % Define time steps to take out of JJAS and keep in MO
        JJASts = (1+(i-1)*120):(105+(i-1)*120);
        MOts = (1+(i-1)*105):(105+(i-1)*105);
        
        % Select correct days
        MOdays(:,:,MOts) = JJASdays(:,:,JJASts);
    end
    
else
    if rem(length(JJASdays(1,1,:)),122) == 0
        % Model is using 365 day calendar
        
        % Find how many years are in time series and make empty array
        yrs = length(JJASdays(1,1,:))/122;
        MOdays = zeros(length(JJASdays(:,1,1)),length(JJASdays(1,:,1)),yrs*107);
        
        for i = 1:yrs
            
            % Define time steps to take out of JJAS and keep in MO
            JJASts = (1+(i-1)*122):(107+(i-1)*122);
            MOts = (1+(i-1)*107):(107+(i-1)*107);
            
            % Select correct days
            MOdays(:,:,MOts) = JJASdays(:,:,JJASts);
        end
        
    else
        disp('WARNING: Calendar of JJAS input data is wrong')
    end
end