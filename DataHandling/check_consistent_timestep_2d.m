function [varout1,timesout1,varout2,timesout2] = check_consistent_timestep_2d(varin1,timesin1,varin2,timesin2)
% [varout1,timesout1,varout2,timesout2] = ...
%     check_consistent_timestep_2d(varin1,timesin1,varin2,timesin2)
%
% This function is designed to ensure data variables from a single model
% are of a consistent length. It was made because one of the files from
% run16 of the UKCP18 GCM ensemble (tasmax) was missing 30 days of data in
% the middle of a time series, so a generic function to correct this was
% written to correct any issues like this that might arise again in the
% future. It does this for 2D variables such as time and yyyymmdd
%
% Outputs:
%   These are the same as the inputs, except the longer variable and set of
%       times steps is reduced to be consistent with the shorter variable
%
% Inputs:
%   varin1 - variable 1 to be processed
%   timesin1 - the 'yyyymmdd' time stamp for variable 1
%   varin2 - variable 2 to be processed
%   timesin2 - the 'yyyymmdd' time stamp for variable 2

%% Find if need to process the data at all
if length(varin1(1,:)) ~= length(varin2(1,:))
    disp('Variables are different length: processing...')
    
    % Convert the time step info into more usable array
    times1 = [];
    times2 = [];
    for i = 1:length(timesin1(1,:))
        times1 = cat(1,times1,str2double(timesin1(:,i)));
    end
    
    for i = 1:length(timesin2(1,:))
        times2 = cat(1,times2,str2double(timesin2(:,i)));
    end
    
    
    %% Find consistent points between longer and shorter timestep
    % If variable 1 is longer
    if length(times1)>length(times2)
        % Shorten variable 1 to match
        [~,correctid1] = intersect(times1,times2,'stable');
        varout1 = varin1(:,correctid1);
        timesout1 = timesin1(:,correctid1);
        
        % Ensure that variable 2 also has the same points
        [~,correctid2] = intersect(times2,times1,'stable');
        varout2 = varin2(:,correctid2);
        timesout2 = timesin2(:,correctid2);
    else
        % If variable 2 is longer
        % Shorten variable 2 to match
        [~,correctid2] = intersect(times2,times1,'stable');
        varout2 = varin2(:,correctid2);
        timesout2 = timesin2(:,correctid2);
        % Ensure that variable 1 also has the same points
        [~,correctid1] = intersect(times1,times2,'stable');
        varout1 = varin1(:,correctid1);
        timesout1 = timesin2(:,correctid1);
    end
    
    
    %% Otherwise the data is already grand and the output is the same as input
else
%     disp('Variables are the same length: grand')
    varout1 = varin1;
    timesout1 = timesin2;
    varout2 = varin2;
    timesout2 = timesin2;
    
end


