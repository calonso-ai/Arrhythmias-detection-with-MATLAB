function features = computeFeatures(arrhythmia, signalArray)

% Compute features from a signal in a segment before an alarm
%
%   features = computeFeatures_new(arrhythmia, signalArray)
%
% Inputs:
%   arrhythmia: the type of arrhythmia the features are computed for
%   signalArray: a struct including the selected signals, sample rate,
%                window for feature computation, beat locations, 
%                signal quality indices, lead indices, and blood pressure 
%                features, if blood pressure signal is present 
%
% Outputs:
%   features: a vector including features for the specified arrhythmia
%   
%   Usage:
%   - signalArray is computed with signalSelection.m
% 
% Written by Linda Eerikäinen
% Last modification by Carlos A. Alonso, May, 2017

if strcmp(arrhythmia, 'Bradycardia')
    features = zeros(1,3);
elseif strcmp(arrhythmia, 'Tachycardia')
    features = zeros(1,2);
elseif strcmp(arrhythmia, 'Ventricular_Tachycardia')
    features = zeros(1,3);

end

Fs = signalArray.fs;
segment_end=Fs*5*60; % alarm position
segment_beginning=segment_end-Fs*signalArray.window+1; % a window before the alarm

% select the beats in the segment
n_ecg_beats=intersect(find(signalArray.rpeaks>=segment_beginning), ...
    find(signalArray.rpeaks<=segment_end));
n_pulse=intersect(find(signalArray.pulse_positions>=segment_beginning), ...
    find(signalArray.pulse_positions<=segment_end));
n_match=intersect(find(signalArray.matched_beats>=segment_beginning), ...
    find(signalArray.matched_beats<=segment_end));

hr_max=NaN(3,1);
max_rr=NaN(3,1);
low_hr_beats = NaN(3,1);
high_hr_beats = NaN(3,1);
low_hr_5 = NaN(3,1);
high_hr_17 = NaN(3,1);

% calculate the heart rate features from the ECG
if length(n_ecg_beats)>=2
    [low_hr_beats(1), low_hr_5(1), high_hr_beats(1), high_hr_17(1), hr_max(1), max_rr(1)] = ...
        hrFeatures(n_ecg_beats, signalArray.rpeaks, segment_end, Fs, 0);
end

% calculate the heart rate features from the pulsatile waveform
if length(n_pulse)>=2
      [low_hr_beats(2), low_hr_5(2), high_hr_beats(2), high_hr_17(2), hr_max(2), max_rr(2)] = ...
    hrFeatures(n_pulse, signalArray.pulse_positions, segment_end, Fs, signalArray.pulse_delay);
    ibis=diff(signalArray.pulse_positions(n_pulse))/Fs;
    ibis=[ibis; ((segment_end-signalArray.pulse_positions(n_pulse(end)) ...
        +signalArray.pulse_delay)/Fs)];
    
else
    ibis = [];
end


% calculate the heart rate features from the matched beats
if length(n_match)>=2
    [low_hr_beats(3), low_hr_5(3), high_hr_beats(3), high_hr_17(3), hr_max(3), max_rr(3)] = ...
        hrFeatures(n_match, signalArray.matched_beats, segment_end, Fs, 0);
end

% In case of ventricular flutter/fibrillation, compute standard deviation
% of inter-beat intervals form pulsatile waveform


% In case of asystole or ventricular flutter/fibrillation, compute features
% from blood pressure
if strcmp(arrhythmia, 'Ventricular_Tachycardia')
    abp_min = NaN;
    if ~isempty(signalArray.abp_features)
        abp_min = nanmin(signalArray.pulsatile(segment_beginning:segment_end));
        if abp_min < 10
            abp_min = NaN;
        end
    end
    
end

% In case of ventricular tachycardia, search beats that differ from average
% QRS complex in the signal
if strcmp(arrhythmia, 'Ventricular_Tachycardia')
    [diffPeaks, checkMatrix] = differentBeatDetection(signalArray.ecg, Fs);
    checkSum = sum(checkMatrix,2);
    diffPeaks = diffPeaks(checkSum >= 2);
    n_diff_peaks =intersect(find(diffPeaks>=segment_beginning),find(diffPeaks<=segment_end));
    num_diff_peaks = length(n_diff_peaks);
end


if strcmp(arrhythmia, 'Bradycardia')
    features(1) = low_hr_5(1);
    features(2) = low_hr_5(2);
    features(3) = low_hr_beats(2);
elseif strcmp(arrhythmia, 'Tachycardia')
    features(1) = high_hr_17(1);
    features(2) = high_hr_beats(1);
elseif strcmp(arrhythmia, 'Ventricular_Tachycardia')
    features(1) = hr_max(1);
    features(2) = hr_max(2);
    features(3) = num_diff_peaks;
end


end


