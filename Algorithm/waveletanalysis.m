function [waveletOutput, wavelet] = waveletanalysis(dataBuffer, startOfSegment, endOfSegment, Fs, Fc, initialize, wavelet)
%WAVELETANALYSIS calculates the multiplication of Mexican hat wavelet
%analyses at two different scales
% 
%   nrOfLevels defines the number of wavelet
%   analyses to be performed on each of the electrode input signals as well
%   as define which lengths the wavelet will have (in case of 3 -> 26, 51,
%   101 samples respectively with an accuracy of 5)
%   
%   - waveletOutput   : array(integer)
%   - lengthOfWavelet : integer
%
%   Written by Michiel Rooijakkers
    
    if initialize
        wavelet = CalculateMexHatWavelet(Fs,Fc);
    end
    
    [waveletOutput] = abs(CalculateWaveletAnalysis(dataBuffer, wavelet, startOfSegment, endOfSegment));
end


function [wavelet] = CalculateMexHatWavelet(Fs,Fc)
%CALCULATEMEXHATWAVELET calculates the Mexican Hat wavelet coefficients 
%based on the given sampling frequency (Fs) and Cut-off frequency (Fc)
    
    baseWidth = (Fs * 0.2516) / Fc;
    limit = round(4*baseWidth) / baseWidth;
    x = -limit:1/baseWidth:limit;  
    
    % Calculate right half of wavelet
    wavelet = ((1-x.^2) .* exp(-(x.^2)/2))';
    % Make total energy of wavelet one
    wavelet = wavelet ./ (ones(size(wavelet)) * sum(abs(wavelet)));
end


function [waveletOutput] = CalculateWaveletAnalysis(dataBuffer, wavelet, startOfSegment, endOfSegment)
%CALCULATEWAVELETOUTPUT calculates the Convolution of the wavelet with the 
%signal and cut off the excess samples at start and end
%
% - waveletOutput     : array(integer)

    % Put data in a single consecutive block
    if startOfSegment < endOfSegment
        dataToAnalyse = dataBuffer(startOfSegment:endOfSegment);
    else
        dataToAnalyse = dataBuffer(startOfSegment:end);
        dataToAnalyse(end+1:end+endOfSegment) = dataBuffer(1:endOfSegment);
    end
        
    % Do convolution
    [tmpWaveletOutput] = conv(dataToAnalyse(1,:), wavelet, 'same');
    % Initialize values
    nrOfSamplesInAnalysis = length(dataToAnalyse);
    offset = (length(wavelet)-1)/2;
    % Put output result into place
    waveletOutput = zeros(1, nrOfSamplesInAnalysis);
    waveletOutput(1+offset:end-offset) = tmpWaveletOutput(1+offset:end-offset);
end