function [HD] = Humidex(T,q,p)
% [W] = Humidex(T,q,p)
% 
% This function calculates the Humidex calculation of Zhao et al. 2015 
% using global and regional data. Vapour pressure (VP is calculated using
% the equation from:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html)
% 

%% Calculate vapour pressure (VP)
VP = (q.*p)./(0.378.*q + 0.622);

%% Calculate Humidex (HD)
HD = T + VP*0.555 - 5.5;

