function [THI] = ThermalHeatIndex(T,q,p)
% [THI] = ThermalHeatIndex(T,q,p)
% 
% Calculate the Thermal Heat Index for agricultural heat stress modelling,
% using the equation in Dunn et al. 2014 (eq. 1), and conversions between
% vapour pressure and relative humidity from:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
% 
%   THI = (1.8*T + 3.2) - (0.55 - 0.0055*RH) * (T - 26.8)
% 
% Outputs:
%   THI = Thermal Heat Index
% 
% Inputs:
%   T = Temperature (�C) - could be max., min. or mean are required
%   q = Specific humidity
%   p = Surface pressure (NOTE: this should be calculated using the same
%       temperaure field (i.e. max., min. or mean)

%% Calculate relative humidity
% First,  Vapor Pressure
VP = VapourPressure(q,p);

% % Then dewpoint temperature
% Td = (243.5 * log (VP/6.112)) ./ (17.67 - log(VP/6.112));

% Then saturation vapour pressure
SVP = 6.112 * exp((17.67*T)./(T+243.5));

% Finally, relative humidity
RH = 100 * VP./SVP;

%% Calculate THI
THI = (1.8*T + 3.2) - (0.55 - 0.0055*RH) .* (T - 26.8);

