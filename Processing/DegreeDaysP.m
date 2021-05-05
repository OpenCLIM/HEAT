function [DDp] = DegreeDaysP(data,p_thresh)
% [DDp] = DegreeDays(data,p_thresh)
% 
% Calculate the number of degree days exceeding a given percentile
% threshold in a dataset. 


% Find the relevant percentile threshold
thresh = prctile(data,p_thresh,3);

% Find days that exceed thresh
days = data > thresh;

% Calculate number of degree days
DDp = sum((data-thresh).*days,3);
