function [HD] = HumidexVP(T,VP)
% [W] = HumidexVP(T,VP)
% 
% This function calculates the Humidex calculation of Zhao et al. 2015 
% using VP data, rather than q and p.
% 

%% Calculate Humidex (HD)
HD = T + VP*0.555 - 5.5;

