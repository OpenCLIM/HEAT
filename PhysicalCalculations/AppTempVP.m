function [AT] = AppTempVP(T,VP)
% [AT] = AppTempVP(T,VP)
% 
% This function calculates the Apparent Temperature calculation of Zhao et al. 2015 
% from VP, rather than q and p.
% 

%% Calculate SWBGT (W)
AT = 0.92.*T + 0.22.*VP - 1.3;

