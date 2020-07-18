function [startOfSegment segmentLength thresholdValue retryCount initialize] = segmentSelection(startOfSegment, segmentLength, thresholdValue, retryCount, peakPosition, ~, Fs, initialize, g)
    if ~isnan(peakPosition)
    % A peak was found in the last interval AND peak not to close to
    % end of segment

        if initialize
            initialize = 0;
            RRInterval = g.RR_MAX;
        else
            RRInterval = peakPosition + g.HALF_BLANK_PERIOD;
        end

        startOfSegment = peakPosition + startOfSegment + g.HALF_BLANK_PERIOD;
        
        if retryCount == g.MAX_TRY-1
            segmentLength = g.RR_MAX;
        else
            segmentLength = round(1.7*RRInterval);
        end
        
        retryCount = 0;
    else
    % No peak was found in the last interval
        thresholdValue = thresholdValue / 2;
        maxSegmentLength   = 2*Fs;

        if retryCount == g.MAX_TRY-1
        % Try somewhere else
            thresholdValue = thresholdValue * (2^retryCount);
            startOfSegment = startOfSegment + Fs;
            segmentLength  = g.RR_MAX;
            retryCount = 1;
        else
        % Retry
            if retryCount == g.MAX_TRY - 2
            % Last Try (2nd retry)
                nrOfSamplesToAdd = (maxSegmentLength - segmentLength);
            else
            % 1st retry
                nrOfSamplesToAdd = round((maxSegmentLength - segmentLength)/2);
            end
            segmentLength = segmentLength + nrOfSamplesToAdd;

            retryCount = retryCount + 1;
        end
    end

    segmentLength = min(max(g.MIN_SEGMENT_LENGTH, segmentLength), g.RR_MAX);
end