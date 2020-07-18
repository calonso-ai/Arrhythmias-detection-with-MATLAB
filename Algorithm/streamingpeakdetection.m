function [out] = streamingpeakdetection(inputData, Fs, HRLimits, PLOT, Fc, blockSize)
%STREAMINGPEAKDETECTION determines the position of all R-peaks in the 1st
%lead of the given inputData in # samples relative to the start of the
%data. 
% The used method is a multiplication of 1 CWT, followed by a peak
% detection with a variable threshold dependent on the local SNR.
%
% Date:   28-3-2011
% Author: Michiel Rooijakkers
% Last modified: Linda Eerikainen, 8-4-2015 
%
% Input arguments:
%  - inputData   : [Array] of input samples
%  - Fs          : [int] sampling frequency of inputData
%  - HRLimits    : [int Array] Upper and lower accepted HR 
%                    (Optional - Default = [36 256])
%  - PLOT        : [bool] Plot results
%                    (Optional - Default = false)
%  - Fc          : [double] indicating high 3dB filter frequency.
%                    (Optional - Default = 18.5)
%  - blockSize   : [int] Size of blocks used in processing the data (in
%                    number of samples.
%                    (Optional - Default = 1024)
%
% Output arguments as part of a structure:
%  - peakPositionArray      : [Array] of found peak positions in nr of
%                               samples since start of file.
%  - waveletOutputArray     : [Array] containing full Wavelet analysis
%                               output result. (Same size as inputData)
%  - thresHoldValueArray    : [Array] Value (amplitude) of all used
%                               thresholds.
%  - startOfSegmentArray    : [Array] Time stamp of start of each segment
%                               in number of samples sice start of file.
%  - endOfSegmentArray      : [Array] Time stamp of end of each segment
%                               in number of samples sice start of file.
    
    if nargin < 6; g.BLOCK_SIZE  = 256; else g.BLOCK_SIZE = blockSize; end;
    if nargin < 5; Fc = 18.5; end;
    if nargin < 4; g.PLOT = false; else g.PLOT = PLOT; end
    if nargin < 3; HRLimits = [36 210]; end
    
    g.NR_OF_SAMPLES        = length(inputData);
    
    if g.PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialization of output arguments
        out.waveletOutputArray           = nan(1,g.NR_OF_SAMPLES);
        out.thresHoldValueArray          = nan(1,0);
        out.startOfSegmentArray          = nan(1,0);
        out.endOfSegmentArray            = nan(1,0);
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.peakPositionArray                = nan(1,0);
    out.SNR                              = nan(1,0);
    out.threshold                        = nan(1,0);
    out.waveMean                         = nan(1,0);
    out.waveMedian                       = nan(1,0);
    out.waveStd                          = nan(1,0);
    
    % Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g.BPM_MIN            = HRLimits(1);
    g.BPM_MAX            = HRLimits(2);
    g.RR_MIN             = floor(Fs * 60/g.BPM_MAX);
    g.RR_MAX             = ceil(Fs * 60/g.BPM_MIN);
    g.DATA_BUFFER_SIZE   = max(g.BLOCK_SIZE*2, 4096);
    g.MIN_SEGMENT_LENGTH = round(Fs/2);
    g.MAX_TRY            = 3;
    g.BASE_RATIO         = 1/3;
    g.HALF_BLANK_PERIOD  = round(g.RR_MIN/2);
    
    dataBuffer           = nan(1, g.DATA_BUFFER_SIZE);

    retryCount           = 0;
    startOfBlock         = 0;
    startOfSegment       = 0;
    endOfSegment         = -1;
    segmentLength        = g.RR_MAX;
    wavelet              = nan;

    ratio                = g.BASE_RATIO;
    minThresholdValue    = 0;

    sampleThreshold      = g.RR_MAX;
    initialize           = 1;
    Fp                   = Fc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    while sampleThreshold + g.BLOCK_SIZE < g.NR_OF_SAMPLES
        % Fill dataBuffer while sampleThreshold has not been exceeded
        while endOfSegment < sampleThreshold
            dataBuffer(mod(startOfBlock, g.DATA_BUFFER_SIZE) + 1 : mod((startOfBlock + g.BLOCK_SIZE - 1), g.DATA_BUFFER_SIZE) + 1 ) ...
                         = inputData(:, startOfBlock+1 : startOfBlock+g.BLOCK_SIZE);
            startOfBlock = startOfBlock + g.BLOCK_SIZE;
            endOfSegment = endOfSegment + g.BLOCK_SIZE;
        end
        
        % Convert segment limits to data buffer indices
        relativeStartOfSegment = mod(startOfSegment, g.DATA_BUFFER_SIZE) + 1;
        relativeEndOfSegment = mod(sampleThreshold, g.DATA_BUFFER_SIZE) + 1;
        
        
        % Perform wavelet analysis
        [waveletOutput wavelet] = waveletanalysis(dataBuffer, relativeStartOfSegment, relativeEndOfSegment, Fs, Fp, initialize, wavelet);
        
        minThresholdValue = (1/8)*(mean(waveletOutput)*2) + (7/8)*minThresholdValue;
        if initialize; thresholdValue = 0.25 * max(waveletOutput,[],2); end
        
        
        [peakPosition thresholdValue] = rpeakdetection(waveletOutput, thresholdValue, ratio, minThresholdValue, g.RR_MIN);
        
        
        if isnan(peakPosition); ratio = g.BASE_RATIO;
        else
            [ratio] = calculateSNR(abs(waveletOutput), peakPosition, Fs);
            wavelet_mean = mean(waveletOutput);
            wavelet_median = median(waveletOutput);
            wavelet_std = std(waveletOutput);
        end
                
        
        waveletLength = length(wavelet);
       if g.PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            out.waveletOutputArray(startOfSegment+round(waveletLength/2):sampleThreshold+1) = waveletOutput(round(waveletLength/2):end);
            out.thresHoldValueArray(end+1)      = thresholdValue;
            out.startOfSegmentArray(end+1)      = startOfSegment;
            out.endOfSegmentArray(end+1)        = sampleThreshold - waveletLength + 1;
        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isnan(peakPosition); out.peakPositionArray(end+1) = startOfSegment + peakPosition; 
            out.SNR(end+1) = ratio;
            out.threshold(end+1) = thresholdValue;
            out.waveMean(end+1) = wavelet_mean;
            out.waveMedian(end+1) = wavelet_median;
            out.waveStd(end+1) = wavelet_std;
            
        end
        
        
        [startOfSegment segmentLength thresholdValue retryCount initialize] ...
            = segmentSelection(startOfSegment, segmentLength, thresholdValue, ...
              retryCount, peakPosition, minThresholdValue, Fs, initialize, g);
        sampleThreshold = startOfSegment + segmentLength;
    end
    
    
% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if g.PLOT
        figure;
        out.waveletArray = nan;
        plotSignals(inputData(1,:),out.peakPositionArray(1,:),out.waveletOutputArray(1,:),out.thresHoldValueArray(1,:),nan,out.startOfSegmentArray(1,:),out.endOfSegmentArray(1,:),Fs, out.waveletArray);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%