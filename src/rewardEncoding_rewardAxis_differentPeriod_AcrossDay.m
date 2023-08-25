close all; clear all;
[dates, rootFolder, axisOutputFileName, processedFileName, analysisBin, params] = pp_getParam("Dataset", "all", "timeperiod", "TO", "neuralDataType", "whole-stitched");
ndates = length(dates);
axisOutputFileName = "../interim/whole-stitched_TO_rewardAxis_all.mat";
outputFolder = rootFolder + "/rewardAxis_vs_parameters_TOrewardAxis_whole/"; makeDir(outputFolder);
outputFolder1DProjection = rootFolder + "/rewardAxisProjection_TOrewardAxis_whole/"; makeDir(outputFolder1DProjection);


trialDataAll = struct.empty;
trialNumBegin = 0;
for d = 1:ndates
    date = dates(d);
    load(processedFileName(d));
    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    
    trialNum = [trialData.trial];
    behavData = behavData(ismember([behavData.trial], trialNum));
    for i = 1:length(trialData)
        trialData(i).reactionTime = behavData(i).reactionTime;
        trialData(i).delayTime = behavData(i).delayTime;
        trialData(i).dayLabel = d;
        trialData(i).newTrial = trialData(i).trial + trialNumBegin;
    end
    trialDataAll = cat(1, trialDataAll, trialData);
    trialNumBegin = trialNumBegin + max(trialNum);
end
trialData = trialDataAll;


% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
reactionTimes = [trialData.reactionTime];
[rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants()
dayLabels = [trialData.dayLabel]; days = unique(dayLabels); ndays = length(days);
trialNums = [trialData.newTrial];
trialStatusLabels = [trialData.trialStatus];
delayTimes = [trialData.delayTime]; delayTimes = round(delayTimes / 50) * 50;


% Get data and avg within dir x rew
% neuralActivity = cat(3, trialData(:).firingRates);
% [nneurons, ntimeBins, ntrials] = size(neuralActivity);
% neuralData = nan(ntrials, nneurons);
% for i = 1:ntrials
%     movementOnset = int16(reactionTimes(i)) + 200;
%     neuralData(i, :) = squeeze(mean(neuralActivity(:, movementOnset+0:movementOnset+200, i), 2));
% end


% Get data and avg within dir x rew
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);

load(axisOutputFileName)


%% Reward Axis
neuralData_onDPC = (neuralData - meanFiringRates) * wR(:,1);

VS_reward_size(neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis")
VS_delay_time(neuralData_onDPC', rewardLabels, difficultyLabels, delayTimes, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis")
VS_direction(neuralData_onDPC', rewardLabels, difficultyLabels, directionLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis")
VS_trial(trialNums, neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis", ...
            "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1650)

%% trajectory
addpath(genpath("./utils/processing/"));
neuralActivity_onDPC = zeros(size(neuralActivity, 2), ntrials);
[neuralActivity_smoothed, ~, ~, ~] = kernelSmooth(double(neuralActivity));
for i = 1:ntrials
    trajectory = squeeze(neuralActivity_smoothed(:, :, i))';
    neuralActivity_onDPC(:, i) = (trajectory - meanFiringRates) * wR(:,1);
end
VS_time(neuralActivity_onDPC, rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis", "Ylim", [-150 50], ...
    "Timeperiod", params.timeperiod, "Xlim", params.Xlim, "DiffLegendPos", "southwest", "RewLegendPos", "northeast");

%% stretching
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);

if params.timeperiod ~= "RF"
    X = zeros(ndays, ndirections); Y = zeros(ndays,ndirections);
    for i = 1:ndays
        for j = 1:ndirections
            curIndsXS = dayLabels == days(i) & directionLabels == directions(j) & difficultyLabels == difficulties(1) & rewardLabels == rewards(1); %Tiny
            curIndsXL = dayLabels == days(i) & directionLabels == directions(j) & difficultyLabels == difficulties(1) & rewardLabels == rewards(3); %Tiny
            curIndsYS = dayLabels == days(i) & directionLabels == directions(j) & difficultyLabels == difficulties(2) & rewardLabels == rewards(1); %Huge
            curIndsYL = dayLabels == days(i) & directionLabels == directions(j) & difficultyLabels == difficulties(2) & rewardLabels == rewards(3); %Huge
            X(i,j) = mean(neuralData_onDPC(curIndsXL)) - mean(neuralData_onDPC(curIndsXS));
            Y(i,j) = mean(neuralData_onDPC(curIndsYL)) - mean(neuralData_onDPC(curIndsYS));
        end
    end

    rankTest2d(X(:), Y(:), "XLabel", "Reward axis range, Tiny targets (Spikes/s)", ...
        "YLabel", "Reward axis range, Huge targets (Spikes/s)", ...
        "OutputFolder", outputFolder+"vs-stretch-SvsL");
else
    X = zeros(1, ndirections); Y = zeros(1,ndirections);
    for j = 1:ndirections
        curIndsXS = directionLabels == directions(j) & difficultyLabels == difficulties(1) & rewardLabels == rewards(1); %Tiny
        curIndsXL = directionLabels == directions(j) & difficultyLabels == difficulties(1) & rewardLabels == rewards(3); %Tiny
        curIndsYS = directionLabels == directions(j) & difficultyLabels == difficulties(2) & rewardLabels == rewards(1); %Huge
        curIndsYL = directionLabels == directions(j) & difficultyLabels == difficulties(2) & rewardLabels == rewards(3); %Huge
        X(j) = mean(neuralData_onDPC(curIndsXS)) - mean(neuralData_onDPC(curIndsXL));
        Y(j) = mean(neuralData_onDPC(curIndsYS)) - mean(neuralData_onDPC(curIndsYL));
    end

    rankTest2d(X(:), Y(:), "XLabel", "Reward axis range, Tiny targets (Spikes/s)", ...
        "YLabel", "Reward axis range, Huge targets (Spikes/s)", ...
        "OutputFolder", outputFolder+"vs-stretch-SvsL");
end


