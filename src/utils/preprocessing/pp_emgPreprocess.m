function [EMGData, EMGMetrics, ECGData] = pp_emgPreprocess(trialData, analogData, outputFormat, plotoutputFolder)
% This function preprocesses the EMG data, including filtering,
% rectification, and downsampling. Then, we get the segment applied to each
% trial in trialData and store it accordingly.
%
% Inputs:
% - analogData: structure with 4 fields:
%   - channelLabels: names for each of the analog channels. For this, we
%   use 'EMG_1' through 'EMG_8'
%   - fs: scalar, sampling frequency of the data
%   - time: [ntime x 1] vector with the time (relative to task CPU) of the
%   analog data
%   - data: [ntime x nchannels] analog data
% - synchInfo: structure with many fields of synch information, but
% for this, all we use is the "taskSynchTrialTimes", which indicates when
% state 1 (center target appearance) first occurred for each trial ([ntrials
% x 1])
% - trialData: [ntrials x 1] structure with the trial parameter data and
% kinematics
%
% Outputs:
% - EMGdata: same as input, but now has 3 more fields:
%   - EMG: [ntime x nmuscles] matrix with EMG values
%   - muscleNames: [nmuscles x 1] cell array with the corresponding 
%     muscle names for each
%   - goodEMGData: [nmuscles x 1] boolean array indicating if the trial's 
%   data for each muscle is good/does not have artifacts or issues
% - EMGMetrics: a structure indicating signal quality for each muscle
%   - baseline = [nmuscles x 1]
%   - maxSignalTuningCurve_mean = [nmuscles x ndirections+1]
%   - maxSignalTuningCurve_std = [nmuscles x ndirections+1]
%   - maxSNR = [nmuscles x 1] (peak avg activity)/(baseline)


%% Filtering
signalData = analogData;
muscleLabel = ["ADel", "LBic", "PDel", "Trap", "Tric"];
fs = 10000;
new_fs = 1000;
zerotime = signalData.time(1); %second

% separate signal by trial
% Step 1 prepare (downsampled length, 5) double array
% Step 2 smoothed and put into the array
% Step 3 separate by trial from synchInfo
% Step 4 load Property from Hand data (reward, direction, success/failure, Handmovement)
EMGs = signalData.data(:, 1:5);
ECG = signalData.data(:, 8);
preprocessedEMGs = zeros(ceil(size(EMGs, 1)/10), 5); % Step 1 prepare (downsampled length, 5) double array
new_fs = 1000;

denoisedECG = pp_ecgPreprocessing(ECG,fs);
for i=(1:length(muscleLabel))
    idx = find(muscleLabel==signalData.EMGMuscleNames(i)); % preventing muscle label misalignment

    baselineRemovedEMG = double(EMGs(:, i));
    for s=(1:8)
        d = designfilt('bandstopiir', 'filterOrder', 2, ...
                        'HalfPowerFrequency1', 60*s-1, 'HalfPowerFrequency2', 60*s+1, ...
                        'DesignMethod', 'butter', 'SampleRate', fs);
        baselineRemovedEMG = filtfilt(d, baselineRemovedEMG);
    end
    baselineRemovedEMG = downsample(baselineRemovedEMG,round(fs/new_fs));
    [ ECGremovedEMG, ~] = pp_ecgRemovalFilter(baselineRemovedEMG, denoisedECG, baselineRemovedEMG(1:100*new_fs), denoisedECG(1:100*new_fs),new_fs);
    bandpassedEMG = bandpass(ECGremovedEMG, [20, 450], new_fs);
    rectifiedEMG = abs(bandpassedEMG);
    preprocessedEMGs(:, idx) = rectifiedEMG; % Step 2 smoothed and put into the array
    
    % ploting
    if ~isempty(outputFormat)
        plotFilteringEffect(downsample(double(EMGs(:, i)),round(fs/new_fs)), baselineRemovedEMG, new_fs, ...
            "Notch Filter "+muscleLabel(idx), plotoutputFolder+muscleLabel(idx)+"-1-baseline-remove", outputFormat);
        plotFilteringEffect(baselineRemovedEMG, ECGremovedEMG, new_fs, ...
            "ECG removal "+muscleLabel(idx), plotoutputFolder+muscleLabel(idx)+"-2-ecg-remove", outputFormat);
        plotFilteringEffect( ECGremovedEMG, bandpassedEMG, new_fs, ...
            "Bandpass 20-450 "+muscleLabel(idx), plotoutputFolder+muscleLabel(idx)+"-3-bandpass-20to450", outputFormat);
        plotFilteringEffect(bandpassedEMG, rectifiedEMG, new_fs, ...
            "Rectify "+muscleLabel(idx), plotoutputFolder+muscleLabel(idx)+"-4-rectify", outputFormat);
    end
end
trialDataEMGRaw = struct.empty(0);
ECGData = struct.empty(0);
for t=(1:length(trialData))
    startTime = ceil((trialData(t).taskSynchTrialTime-zerotime)*1000 + min(trialData(t).time));
    endTime = ceil((trialData(t).taskSynchTrialTime-zerotime)*1000 + max(trialData(t).time));
    trialDataEMGRaw(t).emg = preprocessedEMGs(startTime:endTime, :);
    trialDataEMGRaw(t).time = trialData(t).time;
    trialDataEMGRaw(t).prop.direction = trialData(t).directionLabel;
    trialDataEMGRaw(t).prop.reward = trialData(t).rewardLabel;
    trialDataEMGRaw(t).prop.stateTransition = trialData(t).stateTable;

    ECGData(t).signal = denoisedECG(startTime:endTime, :);
end

%% Normalizing
emgRest = preprocessedEMGs(1:120*new_fs, :);
[EMGData, EMGMetrics] = pp_emgNormalization(trialDataEMGRaw, emgRest, muscleLabel);

%% ploting
dataLength = 0;
for i=(1:length(trialDataEMGRaw))
    stateTransition = trialDataEMGRaw(i).prop.stateTransition;
    if all(ismember([3 4 5 6], stateTransition(1,:))) == 1
        dataLength = dataLength + 1;
    end
end
Y = zeros(dataLength, 5);
s = 0;
for i=(1:length(trialDataEMGRaw))
    stateTransition = trialDataEMGRaw(i).prop.stateTransition;
    if all(ismember([3 4 5 6], stateTransition(1,:))) == 1
        s = s + 1;
        GoCueTime = stateTransition(2, find(stateTransition(1, :)==4));
        EMGaroundGoCue = EMGData(i).signal(GoCueTime-200:GoCueTime+600, :); % -200ms ~ +600ms at GoCue
        Y(s, :) = mean(EMGaroundGoCue, 1);
    end
end
if ~isempty(outputFormat)
    disp('start EMG plot')

    figure
    plot(1:9, EMGMetrics.maxSignalTuningCurve_mean, LineWidth=2.5)
    legend(muscleLabel)
    xticks(1:9)
    xticklabels(["0", "45", "90", "135", "180", "225", "270", "315", "hold"])
    ylabel("Average EMG around peak voltage (a.u.)")
    title("Tuning curve of EMG")
    for e=(1:length(outputFormat))
        saveas(gcf, plotoutputFolder+"0-tuning-curve" + "." + outputFormat(e));
        close all
    end


    for c=(1:length(muscleLabel))
        figure
        plot((1:dataLength), Y(:, c))
        xlabel("Trials")
        ylabel("Mean EMG around GC")
        title("Normalized EMG across trials")
        for e=(1:length(outputFormat))
            saveas(gcf, plotoutputFolder + muscleLabel(c) + "-5-normalized-emg-across-trials" + "." + outputFormat(e));
            close all
        end
    end
end