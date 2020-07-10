function [freq,correct_x] = freqbyarea(inputvar,area,binsize,inputvar2)
% [freq] = freqbyarea(inputvar,area,inputvar2)
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
%   freq = frequency (no. of days) for values of inputvar from -9 to 50 °C. 
%       2D if second inputvar is included.
%   correct_x = the correct values/units for the x bins used (as -9°C has
%       been indexed to x = 1, 0°C is x = 10 etc.)
%
% Inputs:
%   inputvar = the climate variable to be re-shaped
%   area = the absolute grid cell areas at the same resolution as inputvar
%       (use 1 if computing for pre-areally averaged data)
%   binsize = the number of bins to break the data into: 1 = 1°C bins 
%       [default]
%   inputvar2 = optional second climate variable if plotting 

%% Set default binsize and calculate scale factor for bins
if ~exist('binsize','var')
    binsize = 1;
end

scalefac = 1/binsize;


%% Find if calculating in 1D or 2D
% If in 1D
if ~exist('inputvar2','var')
    
    % Generate output file
    freq = zeros(60*scalefac,1);
    correct_x = -9:binsize:51;
    correct_x = correct_x(1:length(freq)); % If using decimal bins, correct_x will be too long
    
    %% Go through each grid cell and repeat climate variable value for each full km2 in its grid cell
    for h = 1:length(inputvar(1,1,:))
        for i = 1:length(inputvar(:,1,1))
            for j = 1:length(inputvar(1,:,1))
                
                % Remove values below -9°C (negligible for summer anyway)
                if inputvar(i,j,h) < -9
                    inputvar(i,j,h) = nan;
                end
                
                % Multiply the inputvar by the scalefac to round to the
                % correct no. of decimal places
                val = floor(inputvar(i,j,h)*scalefac);
                
                % Add current gridcell to correct bin of freq array
                if ~isnan(val)
                    freq(val+(9*scalefac + 1)) = freq(val+(9*scalefac + 1)) + round(area(i,j));
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
