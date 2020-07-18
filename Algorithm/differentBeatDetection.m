function [diffPeaks, checkMatrix] = differentBeatDetection(ECG, Fs)

% Search peaks that are different than regular ECG peaks
%
%   [diffPeaks, checkMatrix] = differentBeatDetection(ECG, Fs)
%
% Inputs:
%   ECG: single lead ECG signal
%   Fs: sampling rate of the ECG signal
%
% Outputs:
%   diffPeak: locations of the different R-peaks in ECG
%   checkMatrix: an array to check which conditions differ sufficiently
%                from the average QRS template
%            Col 1: orientation
%                2: height
%                3: width
%                4: ratio of preceding and succeeding peak interval
%   
% Written by: Linda Eerikäinen, 20 June, 2015
% Thresholds updated 14 August, 2015

diffPeaks = [];
checkMatrix = [];

% Detect peaks from the ECG
out = streamingpeakdetection([ECG; ECG(end-2*Fs+1:end)]', Fs);
peaks = out.peakPositionArray;
peaks = peaks(peaks < 300*Fs);
nPeaks = length(peaks);
if nPeaks > 20
    
    % Form a template to represent a regular QRS complex in the current
    % signal
    [template, tempOrientation] = getTemplate(ECG,peaks,Fs);
    if all(isnan(template))
        return;
    end
    % Compute the width of the template
    templateDiff = diff(template);
    [~,locTemp1] = max(templateDiff);
    [~,locTemp2] = min(templateDiff);
    if tempOrientation == 1
        [templatePeak, peakPosition] = max(template);
        widthTemplate = locTemp2-locTemp1;
    elseif tempOrientation == -1
        [templatePeak, peakPosition] = min(template);
        widthTemplate = locTemp1-locTemp2;
    end
    
    % Peak validation
    templateSize = length(template);
    correlation = NaN(nPeaks,1);
    height = NaN(nPeaks,1);
    width = NaN(nPeaks,1);
    orientation = NaN(nPeaks,1);
    rrRatio = NaN(nPeaks,1);
    delta = round(Fs * 0.010);
    
    % Compare the detected peaks to the template
    if peakPosition > delta && (templateSize-peakPosition) > delta
        for i = 1:nPeaks
            if i > 1 && i < nPeaks
                rrRatio(i) = (peaks(i+1)-peaks(i))/(peaks(i)-peaks(i-1));
            end
            rPeak = ECG(peaks(i)-peakPosition+1:peaks(i)+(templateSize-peakPosition));
            rPeak = rPeak - mean(rPeak);
            orientation(i) = sign(sum(rPeak(peakPosition-delta:peakPosition+delta)));
            rPeakDiff = diff(rPeak);
            [~,loc1] = max(rPeakDiff);
            [~,loc2] = min(rPeakDiff);
            if orientation(i) == 1
                height(i) = max(rPeak);
                width(i) = loc2-loc1;
            elseif orientation(i) == -1
                height(i) = min(rPeak);
                width(i) = loc1-loc2;
            end
            % Compute maximum cross-correlation between template and R-peak
            c_temp = xcorr(rPeak,template,'coeff');
            correlation(i) = max(c_temp);
        end
        
        
        corrThres = median(correlation);
        if corrThres < 0.95
            corrThres = 0.95;
        end
        diffPeaks = peaks(correlation < corrThres);
        diffPeaksIdx = correlation < corrThres;
        nDiffPeaks = length(diffPeaks);
        checkMatrix = zeros(nDiffPeaks,4);
        
        % Set thresholds
        widthThres = round(0.3*widthTemplate);
        heightThres = 0.3*templatePeak;
        ratioThres = 0.7;
        
        % Compare the orientation, height, width, and ratio of preceding
        % and succeeding RR intervals to the template
        checkMatrix(:,1) = orientation(diffPeaksIdx)-tempOrientation;
        if nnz(abs(checkMatrix(:,1)) > 0) > 0
            checkMatrix(:,1) = abs(checkMatrix(:,1)) > 0;
        end
        checkMatrix(:,2) = abs(height(diffPeaksIdx) - templatePeak) > heightThres;
        checkMatrix(:,3) = abs(width(diffPeaksIdx) - widthTemplate) > widthThres;
        checkMatrix(:,4) = rrRatio(diffPeaksIdx) < ratioThres;
        
    end
    
end

end

