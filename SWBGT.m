function [W] = SWBGT(T,q,p)
% [W] = SWBGT(T,q,p)
% 
% This function calculates the simplified WBGT calculation of Zhao et al. 2015 
% using global and regional data. Vapour pressure (VP is calculated using
% the equation from:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html)
% 

%% Calculate vapour pressure (VP)
VP = (q.*p)./(0.378.*q + 0.622);

%% Calculate SWBGT (W)
W = 0.567.*T + 0.393.*VP + 3.94;

