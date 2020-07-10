function [cleanoutput] = remove_RCM_errorpoint(input,threshold)
% [cleanoutput] = remove_RCM_errorpoint(input)
% 
% This function tidies UKCP18 RCM data by removing the corrupted point in
% the Irish Sea which is ~0 for all variables. It converts it to NaN.
% 

%% Set the threshold to 0 if it is not provided
if ~exist('threshold','var')
    threshold = 0;
end

%% Create new file to save output
cleanoutput = input;

%% Remove the dud point
if input(33,52) <= threshold
    cleanoutput(33,52) = NaN;
end