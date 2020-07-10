function [freq] = freqbyarea2(inputvar,area,inputvar2)
% [freq] = freqbyarea2(inputvar,area,inputvar2)
%
% As freqbyarea, but does the calculation by tenth degrees, rather than
% whole degrees
% 
% Convert model data (in lat-long or x-y space) into a simple array for
% plotting e.g. box plots or histograms, where the frequency is based upon
% area. This is a second script/method for doing this and is much more
% efficient.
%
% Also works with two input variables if interested in plotting a bivariate
% histogram (e.g. vapour pressure and temperature on different axes).
%
% Outputs:
%   freq = frequency (no. of km2 days) for values of inputvar from 0:50. 2D
%       if second inputvar is included.
%
% Inputs:
%   inputvar = the climate variable to be re-shaped
%   area = the absolute grid cell areas at the same resolution as inputvar
%   inputvar2 = optional second climate variable if plotting 

if ~exist('inputvar2','var')
    
    freq = zeros(610,1);
    
    %% Go through each grid cell and repeat climate variable value for each full km2 in its grid cell
    for h = 1:length(inputvar(1,1,:))
        for i = 1:length(inputvar(:,1,1))
            for j = 1:length(inputvar(1,:,1))
                
                % Round the climate input at the current grid cell to the
                % nearest 10th, then add 10 to allow for negative values,
                % then multiply by 10 to convert to whole numbers for
                % indexing the output... Simples.
                val = round(inputvar(i,j,h)+10,1)*10;
                
                % If a value is less than -10, make it -10
                if val < 1
                    val = 1;
                else
                    if val > 610
                        val = 610;
                    end
                end
                
                
                if ~isnan(val)
                    freq(val) = freq(val) + area(i,j);
                end
            end
        end
    end
    
else
    
    freq = zeros(51,51);
    
    %% Go through each grid cell and repeat climate variable value for each full km2 in its grid cell
    for h = 1:length(inputvar(1,1,:))
        for i = 1:length(inputvar(:,1,1))
            for j = 1:length(inputvar(1,:,1))
                
                val = floor(inputvar(i,j,h));
                val2 = floor(inputvar2(i,j,h));
                
                
                if val < 1
                    val = 1;
                end
                if val2 < 1
                    val2 = 1;
                end
                
                
                if ~isnan(val)
                    freq(val,val2) = freq(val,val2) + round(area(i,j));
                end
            end
        end
    end
end
