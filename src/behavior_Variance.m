%% this function is to calculate signal variance and noise variance
% of behaviors

% get the reaction time and peak speed across trials
% (1) signal variance
% get the mean reaction time and peak speed in each reward and 
% difficulty condition
% calculate the difference between Small reward and Large reward
% (2) noise variance
% get the variance of reaction time and peak speed in each reward and
% difficulty condition

clear all; close all;
% input parameters
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
ndates = length(dates);
Yall_SV = cell(ndates, 1);
Yall_NV = cell(ndates, 1);

for d = 1:ndates
    date = dates(d);
    rootDir = "../";
    load(rootDir+"data/processed/non-stitched/2022"+date+"_GC_200_600_3456.mat");
    outputFolder = "../results/202306w3-results-summary-all/behavior/behvaior_variance/"; makeDir(outputFolder);

    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    trialNum = [trialData.trial];

    behavData = behavData(ismember([behavData.trial], trialNum));
    for i = 1:length(trialData)
        trialData(i).reactionTime = behavData(i).reactionTime;
        trialData(i).peakSpeed = behavData(i).peakSpeed;
        trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
    end

    % Get labels
    directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
    rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
    difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).'; % Huge=1, Tiny=0
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
    diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    reactionTimes = [trialData.reactionTime];
    peakSpeeds = [trialData.peakSpeed];
    delayTimes = [trialData.delayTime];

    % get the mean reaction time and peak speed in each reward and
    meanReactionTime = zeros(ndifficulties, nrewards, ndirections);
    meanPeakSpeed = zeros(ndifficulties, nrewards, ndirections);
    varReactionTime = zeros(ndifficulties, nrewards, ndirections);
    varPeakSpeed = zeros(ndifficulties, nrewards, ndirections);

    for r = 1:nrewards
        for j = 1:ndifficulties
            for i = 1:ndirections
                curInds = rewardLabels==rewards(r) & difficultyLabels==difficulties(j) & directionLabels==directions(i);
                % Attention: there are NAN values in reaction time and peak speed
                meanReactionTime(j, r, i) = nanmean(reactionTimes(curInds));
                meanPeakSpeed(j, r, i) = nanmean(peakSpeeds(curInds));
                varReactionTime(j, r, i) = nanvar(reactionTimes(curInds));
                varPeakSpeed(j, r, i) = nanvar(peakSpeeds(curInds));
            end
        end
    end

    % calculate signal variance
    signalVarReactionTime = zeros(ndifficulties, ndirections);
    signalVarPeakSpeed = zeros(ndifficulties, ndirections);
    for j = 1:ndifficulties
        for i = 1:ndirections
            signalVarReactionTime(j, i) = meanReactionTime(j, 3, i) - meanReactionTime(j, 1, i);
            signalVarPeakSpeed(j, i) = meanPeakSpeed(j, 3, i) - meanPeakSpeed(j, 1, i);
        end
    end

    % calculate noise variance difference
    noiseVarReactionTime = zeros(nrewards, ndirections);
    noiseVarPeakSpeed = zeros(nrewards, ndirections);
    for r = 1:nrewards
        for i = 1:ndirections
            noiseVarReactionTime(r, i) = varReactionTime(2, r, i) - varReactionTime(1, r, i); % Huge - Tiny
            noiseVarPeakSpeed(r, i) = varPeakSpeed(2, r, i) - varPeakSpeed(1, r, i);
        end
    end

    signalVarances = zeros(2, ndifficulties, ndirections);
    signalVarances(1, :, :) = signalVarReactionTime;
    signalVarances(2, :, :) = signalVarPeakSpeed;
    noiseVarances = zeros(2, nrewards, ndirections);
    noiseVarances(1, :, :) = noiseVarReactionTime;
    noiseVarances(2, :, :) = noiseVarPeakSpeed;

    Yall_SV{d} = signalVarances;
    Yall_NV{d} = noiseVarances;
end


%% analysis of signal variance

Y = zeros(ndates, 2, ndifficulties, ndirections);
for d = 1:ndates
    Y(d, :, :, :) = Yall_SV{d};
end

statusLabel = ["RT", "PS"];
unitLabel = ["(ms)", "(m/s)"];
for status = 1:2
    x = zeros(ndates, ndirections); y = zeros(ndates, ndirections);
    for d = 1:ndates
        for j = 1:ndirections
            x(d,j) = Y(d, status, 1, j);
            y(d,j) = Y(d, status, 2, j);
        end
    end
    rankTest2d(x(:), y(:), "XLabel", statusLabel(status) + " range, Tiny targets " + unitLabel(status), ...
        "YLabel", statusLabel(status) + " range, Huge targets " + unitLabel(status), ...
        "AxisColor","diffColor", "OutputFolder",outputFolder+"behavior-stretch-"+statusLabel(status))
    
    corrHistogram(y(:) - x(:) , "Label", "Huge range - Tiny range", ...
            "OutputFolder", outputFolder+"behavior-stretch-hist"+statusLabel(status), ...
            "Xlim","none");
end

%% analysis of noise variance

Y = zeros(ndates, 2, nrewards, ndirections);
for d = 1:ndates
    Y(d, :, :, :) = Yall_NV{d};
end

statusLabel = ["RT", "PS"];
for status = 1:2
    rewardNames = ["Small", "Medium", "Large"];
    for r = 1:nrewards
        diffxy = squeeze(Y(:, status, r, :));
        corrHistogram(diffxy(:), "Label", "Huge varinace - Tiny varinace:" + rewardNames(r), ...
            "OutputFolder", outputFolder+"behavior-variance-diff-hist-"+statusLabel(status)+"-"+rewardNames(r), ...
            "Xlim","none");
    end
end

close all

