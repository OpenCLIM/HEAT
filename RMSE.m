function [val] = RMSE(var1,var2,area)
% [val] = RMSE(var1,var2,area)
% 
% Calculate the root mean square error between two spatial datasets over a
% given area. This is used during the UKCP18 evaluation with ERA5 and
% HadUK-Grid.
% 
% Outputs:
%   val = the value of RMSE for the entire spatial domain
% 
% Inputs:
%   var1 = the 'correct' data, i.e. WFDEI reanalysis data
%   var2 = the data to compare, i.e. the simulated UKCP data
%   area = the area over which to calculate the RMSE (fractional)
% 

% Calculate the anomaly/error
E = var2-var1;

% Square the error
SE = E.^2;

% Find mean of squared error
MSE = nansum(nansum(SE.*area));

% Calculate the RMSE
val = sqrt(MSE);