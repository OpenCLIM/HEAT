function [VP] = VapourPressure(q,p)
% [W] = Humidex(q,p)
% 
% This function calculates the Vapour Pressure calculated using the equation from:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
% 

%% Calculate vapour pressure (VP)
VP = (q.*p)./(0.378.*q + 0.622);

