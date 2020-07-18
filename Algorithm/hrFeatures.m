function [low_hr_beats, low_hr_5, high_hr_beats, high_hr_17, hr_max, max_rr] = ...
    hrFeatures(n_beats, ann, alarm_time, Fs, delay)

% Compute heart rate features from a signal in a segment before an alarm
%
%   [low_hr_beats, low_hr_5, high_hr_beats, high_hr_17, hr_max, max_rr] = ...
%   hr_features(n_beats, ann, time_alarm, Fs, delay)
%
% Inputs:
%   n_beats: beat indices in  16-second segment before the alarm
%   ann: beat annotations in a signal
%   alarm_time: time of the alarm
%   delay: delay between beats detected from ECG
%   Fs: sample frequency of the signal
%
% Outputs:
%   low_hr_beats: number of consecutive beats with heart rate lower than 40
%                 bpm
%   low_hr_5: minimum heart rate of 5 consecutive beats
%   high_hr_beats: number of consecutive beats with heart rate higher than
%           140 bpm
%   high_hr_17: maximum heart rate of 17 consecutive beats
%   hr_max: maximum heart rate
%   max_rr: maximum RR interval
%
% 
% Written by: Linda Eerikäinen, 1 April, 2015
% Modified tolerance 12 August, 2015, by Linda Eerikäinen
% Last modification by Carlos A. Alonso, May, 2017 

tolerance = 0;

low_hr_min = zeros(length(n_beats)-1,1);
    high_hr_max = zeros(length(n_beats)-1,1);
    for n = 1:length(n_beats)-1
        low_hr = zeros(length(n_beats)-n,1);
        high_hr = zeros(length(n_beats)-n,1);
        % Calculate heart rate from i consecutive beats
        for i=1:length(n_beats)-n
            low_hr(i)=60*Fs/((ann(n_beats(i+n)) ...
                -ann(n_beats(i)))/n);
            high_hr(i)=60*Fs/((ann(n_beats(i+n)) ...
                -ann(n_beats(i)))/n);
        end
        low_hr_min(n) = min(low_hr);
        high_hr_max(n) = max(high_hr);
    end
    % Maximum hear rate in the segment
    hr_max=60*Fs/min(diff(ann(n_beats)));
    max_rr_beats=max(diff(ann(n_beats)))/Fs;
    % Maximum RR interval in the segment
    max_rr=max(max_rr_beats, ((alarm_time-ann(n_beats(end))+delay)/Fs));
    % Search beats with heart rate lower than 40 bpm + tolerance
    low_beats = find(low_hr_min <= 40 + tolerance);
    if ~isempty(low_beats)
        low_hr_beats = max(low_beats)+1;
    else
        low_hr_beats = 0;
    end
    % Search beats with heart rate higher than 140 bpm - tolerance
    high_beats = find(high_hr_max >= 140 - tolerance);
    if ~isempty(high_beats)
        high_hr_beats = max(high_beats)+1;
    else
        high_hr_beats = 0;
    end
    % Heart rate of 17 consecutive beats at maximum and 5 consecutive beats
    % at minimum
    if length(high_hr_max) >= 16
        high_hr_17 = high_hr_max(16);
        low_hr_5 = low_hr_min(4);
    elseif length(low_hr_min) >= 4
        low_hr_5 = low_hr_min(4);
        high_hr_17 = NaN;
    else
        low_hr_5 = NaN;
        high_hr_17 = NaN;
    end
        
end