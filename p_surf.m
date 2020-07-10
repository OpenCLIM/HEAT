function [p_s] = p_surf(psl,T,h)
% [p_s] = p_surf(psl,T,h)
% 
% This function calculates the surface pressure (p_s, mb) at ground level 
% of elevation (h, m), from model output of mean temperature (T, K) and 
% mean sea level pressure (psl, mb) using the hypsometric equation:
% 
%   p_s = psl/exp((g*h)/(R*T))
%
% Where:
% g = acceleration due to gravity (9.80665 m/s^2)
% R = specific gas constant for dry air (287.058 J/kg/K)
% 
% Note: T input in °C will be automatically corrected to K


%% Correct if T is in °C
if T < 200
    T = T+273.15;
end

%% Calculate surface pressure
p_s = psl./exp((9.80665.*h)./(287.058.*T));