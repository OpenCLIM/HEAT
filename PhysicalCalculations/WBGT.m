function [W] = WBGT(T,q,p)
% Calculate the WBGT as used in Andrews et al. 2018
% 

% Calculate Dew Point Temp. from Vapour Pressure
VP = VapourPressure(q,p);
DPT = 243*log(VP/0.6112)/17.67 - log(VP/0.6112);

% Calculate gamma and theta
thet = (4098 * VP)./(DPT + 237.3).^2;
gam = 0.00066*p;

% Calculate the psychometric wet bulb temperature
Tpwb = (gam.*T + thet.*DPT)./(gam + thet);

% Calculate the Wet Bulb Globe Temperature
W = 0.67*Tpwb + 0.33*T;