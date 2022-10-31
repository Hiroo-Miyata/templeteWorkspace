function [normalizedTrialData, EMGMetrics] = emgNormalization(preprocessedTrialData, emgRest, muscleLabel)
%% this function is normalize the emg data 
% normalize method
% 1. get EMG signal of successes
% 2. get averages of peak values of EMG during reach period at each
% direction and that of mean values of EMG during delay period(-200 ~ 0 ms) 
% 3. culculate z-scoring parameter of the 9 datapoints(tuning curve)
% 4. normalize all EMG signal by the parameters 

nmuscles = 5;
% prepare variables for visualization
tuningCurve = zeros(9, nmuscles);
tuningCurveStd = zeros(9, nmuscles);
normalizedParams = zeros(2, nmuscles);


%% get EMG signal of successes 
dataLength = 0;
for i=(1:length(preprocessedTrialData))
    stateTransition = preprocessedTrialData(i).prop.stateTransition;
    if all(ismember([3 4 5 6], stateTransition(1,:))) == 1
        dataLength = dataLength + 1;
    end
end

s = 0;
EMG = zeros(801, nmuscles, dataLength);
directionArray = zeros(1, dataLength);
rewardArray = zeros(1, dataLength);
for i=(1:length(preprocessedTrialData))
    stateTransition = preprocessedTrialData(i).prop.stateTransition;
    if all(ismember([3 4 5 6], stateTransition(1,:))) == 1
        s = s+1;
        GoCueTime = stateTransition(2, find(stateTransition(1, :)==4));
        EMGaroundGoCue = preprocessedTrialData(i).emg(GoCueTime-200:GoCueTime+600, :); % -200ms ~ +600ms at GoCue
        EMG(:,:, s) = EMGaroundGoCue;
        directionArray(s) = preprocessedTrialData(i).prop.direction;
        rewardArray(s) = preprocessedTrialData(i).prop.reward;
    end
end

%% get averages of peak values
for direction=(1:8)
    oneDirectionEMG = EMG(:,:,directionArray==direction);
    movmeanDirectionEMG = movmean(oneDirectionEMG, 100, 1);
    MaxIntensitysAtOneDirection = max(movmeanDirectionEMG, [], 1);
    tuningCurve(direction, :) = mean(MaxIntensitysAtOneDirection, 3);
    tuningCurveStd(direction, :) = std(MaxIntensitysAtOneDirection, 0, 3);
end
tuningCurve(9, :) = mean(EMG(1:200, :, :), [1 3]); % mean at delay period
tuningCurveStd(9, :) = std(EMG(1:200, :, :), 0, [1 3]);
baseline = mean(emgRest, 1);

%% culculate z-index parameters
for channel = (1:nmuscles)
    Y = tuningCurve(:, channel);
    Ymean = mean(Y);
    Ystd = std(Y);
    normalizedParams(:, channel) = [Ymean Ystd];
end

%% normalize all EMG signal by the parameters 
normalizedTrialData = struct.empty(0);
for i=(1:length(preprocessedTrialData))
    normalizedTrialData(i).signal = zeros(size(preprocessedTrialData(i).emg));
    for channel = (1:nmuscles)
        normalizedTrialData(i).signal(:, channel) = (preprocessedTrialData(i).emg(:, channel) - normalizedParams(1, channel)) ./ normalizedParams(2, channel);
    end
    normalizedTrialData(i).goodEMGData = true(1, nmuscles);
end


EMGMetrics = struct();
EMGMetrics.muscleNames = muscleLabel;
EMGMetrics.baseline = baseline;
EMGMetrics.maxSignalTuningCurve_mean = tuningCurve;
EMGMetrics.maxSignalTuningCurve_std = tuningCurveStd;
EMGMetrics.maxSNR = max(tuningCurve, [], 1) ./ baseline;
EMGMetrics.normalizedParams = normalizedParams;