function [W] = SWBGTVP(T,VP)
% [W] = SWBGTVP(T,VP)
% 
% This function calculates the simplified WBGT calculation of Zhao et al. 2015 
% using VP data, rather than q and p.
% 

%% Calculate SWBGT (W)
W = 0.567.*T + 0.393.*VP + 3.94;

