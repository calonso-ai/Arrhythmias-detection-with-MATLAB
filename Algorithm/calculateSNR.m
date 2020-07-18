function [ratio] = calculateSNR(waveletOutput, peakPosition, Fs)
%CALCULATE a number indicating the noise level in the signal around the
%last found peak position.
%
%   Written by Michiel Rooijakkers

%#ok<*UNRCH>
%#ok<*AGROW>
    DEBUG = false;
    
    DELTA       = round(0.075*Fs);
    isoWidth    = round(0.175*Fs);%0.1
    spaceP      = round(0.075*Fs);%0.15
    spaceT      = round(0.075*Fs);%0.15

    isoPS  = 1; 
    isoPE  = peakPosition - spaceP;
    isoTS  = peakPosition + spaceT;
    isoTE  = min(peakPosition+isoWidth+spaceT-1,size(waveletOutput,2));
    peakRS = peakPosition - DELTA;
    peakRP = peakPosition;
    peakRE = min(peakPosition + DELTA,length(waveletOutput));
    
    isoInterval = zeros(1,0);

    % If there is enough space around R-peak, select interval around R-peak
    % and see if enough iso-electrical period is available for Noise level
    % estimation
    if peakRP > 0
        peakInterval = waveletOutput(peakRS:peakRE);

        if isoPS < isoPE
            isoInterval = [isoInterval waveletOutput(isoPS:isoPE)];
        end

        if isoTS < isoTE
            isoInterval = [isoInterval waveletOutput(isoTS:isoTE)];
        end
    end

    SNR = log2(max(peakInterval)/max(isoInterval)); 
    ratio = (6 - SNR) / 8;
    
    % DEBUG PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DEBUG 
        plot(waveletOutput);
        hold on;
        stem(iso,repmat(2*max(waveletOutput),1,4),'r');
        stem(peak,repmat(2*max(waveletOutput),1,3),'g');
        ylim([0 max(waveletOutput)*1.2]);
        ind = (ratio - 1)/5 * 1.2*max(waveletOutput);
        plot([0 400], [ind ind], 'm');
        hold off;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

