function [VP] = DPT2VP(DPT)
% [VP] = DPT2VP(DPT)
% 
% Convert dew point temperature (DPT) to vapour pressure (VP) using the
% equation from: https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html

VP = 6.112*exp((17.67*DPT)./(DPT + 243.5)); 
