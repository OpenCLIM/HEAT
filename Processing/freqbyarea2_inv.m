function [freq_inv] = freqbyarea2_inv(freq)
% [freq_inv] = freqbyarea2_inv(freq)
% 
% Invert the output of freqbyarea2 for use with a Q-Q plot
% 



for i = 1:length(freq)
    if i == 1
        freq_inv = i * ones(round(freq(i)*10),1);
    else
        freq_inv = cat(1,freq_inv,i * ones(round(freq(i)*10),1));
    end
end