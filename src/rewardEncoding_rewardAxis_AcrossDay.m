close all; clear all;
[dates, rootFolder, axisOutputFileName, processedFileName, analysisBin] = pp_getParam("Dataset", "all", "timeperiod", "TO", "neuralDataType", "whole-stitched");
ndates = length(dates);
outputFolder = rootFolder + "/rewardAxis_vs_parameters_wholestitched/"; makeDir(outputFolder);
outputFolder1DProjection = rootFolder + "/rewardAxisProjection_wholestitched/"; makeDir(outputFolder1DProjection);
% analysisBin = (100:200); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(400:600), SU=(200:400)

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
[rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants()
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
trialStatusLabels = [trialData.trialStatus];
delayTimes = [trialData.delayTime]; delayTimes = round(delayTimes / 50) * 50;

% Get data and avg within dir x rew
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);
neuralData_byParameters_mean = nan(ndirections,nrewards,nneurons);
for l = 1:ndirections
    for r = 1:nrewards
            curInds = directionLabels==directions(l) & rewardLabels==rewards(r);
            neuralData_byParameters_mean(l,r,j,:) = mean(neuralData(curInds,:));
    end; clear r
end; clear l

% reward axis of both
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1]));
% meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zR(nrewards, 1) < zR(1, 1)
    wR(:, 1) = -wR(:, 1);
end
meanFiringRates = mean(neuralData, 1);

save(axisOutputFileName, "wR", "eigVls_R", "meanFiringRates");

neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
%% Reward Axis
VS_reward_size(neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis", "Ylim", [-25 20])
VS_delay_time(neuralData_onDPC', rewardLabels, difficultyLabels, delayTimes, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis")
VS_direction(neuralData_onDPC', rewardLabels, difficultyLabels, directionLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis")
VS_trial(trialNums, neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis", ...
            "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1650)