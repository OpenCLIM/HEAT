function [DDa] = DegreeDays(data,abs_thresh)
% [DDa] = DegreeDays(data,abs_thresh)
% 
% Calculate the number of degree days exceeding a given absolute
% threshold in a dataset. 

% Find days that exceed thresh
days = data > abs_thresh;

% Calculate number of degree days
DDa = sum((data-abs_thresh).*days,3);
