function [AT] = AppTemp(T,q,p)
% [AT] = AppTemp(T,q,p)
% 
% This function calculates the Apparent Temperature calculation of Zhao et al. 2015 
% using global and regional data. Vapour pressure (VP is calculated using
% the equation from:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html)
% 

%% Calculate vapour pressure (VP)
VP = (q.*p)./(0.378.*q + 0.622);

%% Calculate SWBGT (W)
AT = 0.92.*T + 0.22.*VP - 1.3;

