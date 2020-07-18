function signalArray = signalSelection(recName)

% Select one ECG lead and one pulsatile waveform from a record based on
% signal quality measures
%
%   signalArray = signalSelection(recName)
%
% Inputs:
%   recName: name of the signal file
%
% Outputs:
%   signalArray: a struct containing
%                   - fs: sample frequency of the signals 
%                   - window: window in which SQIs are computed
%                   - ecg: selected ECG signal
%                   - sqi_features: signal quality features for ECG
%                   - ecg_lead: index for the selected ECG lead
%                   - rpeaks: peak locations in the selected ECG signal
%                   - pulsatile: selected pulsatile waveform signal
%                   - pulse_sqi: signal quality index
%                   - pulse_positions: positions of the pulses detected in
%                     pulsatile signal
%                   - pulse_delay: delay between selectile pulsatile signal
%                     and ECG
%                   - abp_features: arterial blood pressure features 
%                   - matched_beats: locations of beats that match in
%                     several signals
%   Usage:
%   - sqi_features  are obtained using ecgSqiFeatures.m
%   - abp_features  are obtained using abpfeature.m

% 
% Written by: Linda Eerikainen, 3 August, 2015
% Last modification by: Carlos A. Alonso, May, 2017

% Read the signal from a .mat file
[~, signal, Fs, siginfo] = rdmat(recName);

%Resample the signals
Fs=Fs(1);
if Fs~=125
    signal=resample(signal,125,Fs);
    Fs=125;
end

signalArray.fs = Fs;
signalArray.window = 16;

%Replace missing ECG with zeros
signal(isnan(signal(:,1)),1) = 0;
signal(isnan(signal(:,2)),2) = 0;

% Compute signal quality features
[sqi(1).features, peakArray(1).positions, peakArray(1).matches] = ...
    ecgSqiFeatures(signal(:,1), signal(:,2), Fs, signalArray.window);
[sqi(2).features, peakArray(2).positions, peakArray(2).matches] = ... 
    ecgSqiFeatures(signal(:,2), signal(:,1), Fs, signalArray.window);

ppg_idx = [];
abp_idx = [];
mean_sqi_abp = -1;
mean_sqi_ppg = -1;
for i = 3:size(siginfo,2)
    if strfind(siginfo(i).Description, 'PLETH')
        ppg_idx = i;
    elseif strfind(siginfo(i).Description, 'ABP')
        abp_idx = i;
    end
end

% If blood pressure is present, find beats, and compute features and quality
% index
if ~isempty(abp_idx)
    abp = signal(:,abp_idx);
    abp_pulses = wabpC(abp,0,1);
    features_abp = abpfeature(abp,abp_pulses);
    BeatQ = jSQI(features_abp, abp_pulses, abp);
    if ~isempty(BeatQ)
        abp_pulses_in_window = find(abp_pulses > Fs*(300-signalArray.window)+1 & abp_pulses < Fs*300);
        abp_sqi_by_beats=1-sum(BeatQ(abp_pulses_in_window(1:end-1),2:end),2) ...
            /size(BeatQ(abp_pulses_in_window(1:end-1),2:end),2);
        mean_sqi_abp = mean(abp_sqi_by_beats);
    else
        mean_sqi_abp = 0;
    end
end

% If photoplethysmogram is present, find beats, and compute quality index
if ~isempty(ppg_idx)
    ppg = signal(:,ppg_idx);
    y=quantile(ppg,[0.05,0.5,0.95]);%?quantile? algorithm was applied to partition the signal into three quantiles, (0.05, 0.5, 0.95)
    ppg_pulses=wabpC(ppg,0,(y(3)-y(1))/120);
    if ~isempty(ppg_pulses)
        psqi=ppgSQI(ppg,ppg_pulses);
        ppg_pulses_in_window = find(ppg_pulses > Fs*(300-signalArray.window)+1 & ppg_pulses < Fs*300);
        if ~isempty(ppg_pulses_in_window) && length(ppg_pulses)-1 == length(psqi)
            mean_sqi_ppg = mean(psqi(ppg_pulses_in_window(1:end-1)));
        else
            mean_sqi_ppg = 0;
        end
    else
        mean_sqi_ppg = 0;
    end
end


delay_abp = 33; % average ABP delay from ECG in samples
delay_ppg = 61; % average PPG delay from ECG in samples

% Find matching beats between two ECG leads
matches1 = zeros(length(peakArray(1).matches),3);
matches2 = zeros(length(peakArray(2).matches),3);
matches1(:,1) = peakArray(1).matches;
matches2(:,1) = peakArray(2).matches;

% Find matches between ECG and ABP
if ~isempty(abp_idx)
    matches1(:,2) = beatMatch(peakArray(1).positions, abp_pulses-delay_abp, Fs);
    matches2(:,2) = beatMatch(peakArray(2).positions, abp_pulses-delay_abp, Fs);
end

% Find matches between ECG and PPG
if ~isempty(ppg_idx)
    matches1(:,3) = beatMatch(peakArray(1).positions, ppg_pulses-delay_ppg, Fs);
    matches2(:,3) = beatMatch(peakArray(2).positions, ppg_pulses-delay_ppg, Fs);
end

% Sum and sort the found matches
match_sum1 = sum(matches1,2);
match_sum2 = sum(matches2,2);
matches_tot = [peakArray(1).positions(match_sum1 > 0), ...
    peakArray(2).positions(match_sum2 > 0)];
matches_sorted = sort(matches_tot);

%Remove the matching beats that are close (were considereded as match)
matched_beats = matches_sorted(diff(matches_sorted) > 60*Fs/300);


% Select ECG channel and save signal, quality features, lead label, and
% R-peak positions to an array
if isnan(sqi(2).features(1)) %Aqui se como se accede al parametro 1 del vector proporcionado por ecgSqifeatures
    signalArray.ecg = signal(:,1);%que representa el "median local noise of R-peaks"
    signalArray.sqi_features = sqi(1).features;%sqi(1) y sqi(2) son los vectores sqi´s del canal ECG 1 y 2 respectivamente
    signalArray.ecg_lead = 1;
    signalArray.rpeaks = peakArray(1).positions;
elseif sqi(1).features(1) <= sqi(2).features(1) && sqi(1).features(2) > 0
    signalArray.ecg = signal(:,1);
    signalArray.sqi_features = sqi(1).features;
    signalArray.ecg_lead = 1;
    signalArray.rpeaks = peakArray(1).positions;
elseif sqi(2).features(1) < sqi(1).features(1) && sqi(2).features(2) > 0
    signalArray.ecg = signal(:,2);
    signalArray.sqi_features = sqi(2).features;
    signalArray.ecg_lead = 2;
    signalArray.rpeaks = peakArray(2).positions;
elseif sqi(1).features(2) > sqi(2).features(2)
    signalArray.ecg = signal(:,1);
    signalArray.sqi_features = sqi(1).features;
    signalArray.ecg_lead = 1;
    signalArray.rpeaks = peakArray(1).positions;
else
    signalArray.ecg = signal(:,2);
    signalArray.sqi_features = sqi(2).features;
    signalArray.ecg_lead = 2;
    signalArray.rpeaks = peakArray(2).positions;
end

% Select pulsatile waveform and save signal, quality index, lead label,
% pulse positions, and ABP features to an array
if mean_sqi_abp > 0.85
    signalArray.pulsatile = signal(:,abp_idx);
    signalArray.pulse_sqi = mean_sqi_abp;
    signalArray.pulse_lead = abp_idx;
    signalArray.pulse_positions = abp_pulses;
    signalArray.pulse_delay = delay_abp;
    signalArray.abp_features = features_abp;
elseif mean_sqi_abp >= mean_sqi_ppg
    signalArray.pulsatile = signal(:,abp_idx);
    signalArray.pulse_sqi = mean_sqi_abp;
    signalArray.pulse_lead = abp_idx;
    signalArray.pulse_positions = abp_pulses;
    signalArray.pulse_delay = delay_abp;
    signalArray.abp_features = features_abp;
else
    signalArray.pulsatile = signal(:,ppg_idx);
    signalArray.pulse_sqi = mean_sqi_ppg;
    signalArray.pulse_lead = ppg_idx;
    signalArray.pulse_positions = ppg_pulses;
    signalArray.pulse_delay = delay_ppg;
    signalArray.abp_features = [];
end

% Save matched beats to signal array
signalArray.matched_beats = matched_beats;

end