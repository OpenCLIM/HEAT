function [output] = HadUK2UKCP18(input)
% [output] = HadUK2UKCP18(input)
% 
% The HadUK-Grid (12km) dataset and the UKCP18 RCM (12km) datasets do not
% line up exactly. This script sorrects for this, turning the HadUK-Grid
% data (83 x 110 grid cells) to UKCP18 RCM grid cells (82 x 112).

if ndims(input) == 2
    % Create empty array of correct size to fill
    output = zeros(82,112);
    
    % Move most of the grid to the correct place
    output(3:82,2:111) = input(1:80,:);
    
    % Fill edges with empty values (= middle of Atlantic)
    output(1:2,:) = output(3,2);
    output(:,1) = output(3,2);
    output(:,112) = output(3,2);
    
else
    if ndims(input) == 3
        % Create empty array of correct size to fill
        output = zeros(82,112,length(input(1,1,:)));
        % Move most of the grid to the correct place
        output(3:82,2:111,:) = input(1:80,:,:);
        
        % Fill edges with empty values (= middle of Atlantic)
        output(1:2,:,:) = output(3,2,1);
        output(:,1,:) = output(3,2,1);
        output(:,112,:) = output(3,2,1);

    end
end

