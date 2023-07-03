close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"]; % "0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"
fileName = "all";
outputFolderSR = "../results/202306w3-results-summary-"+ fileName + "/behavior_success-rate/";
outputFolderFailureMode = "../results/202306w3-results-summary-"+ fileName + "/behavior_failure-modes/";
outputFolderVigor = "../results/202306w3-results-summary-"+ fileName + "/behavior_vigors/";

%% conbine across days
trialDataAll = struct.empty;
behavDataAll = struct.empty;
trialNumBegin = 0;
for d = 1:length(dates)
    date = dates(d);
    dataFolder = "../data/preprocessed/2022"+date+"/";
    taskfile = dir(dataFolder+ "*" +date+ "*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    taskfile = dir(dataFolder+ "*" +date+ "*_TSK_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);

    trialNum = [trialData.trial];
    trialData = trialData(trialNum > 24);
    trialNum = [trialData.trial];
    behavData = behavData(ismember([behavData.trial], trialNum));
    
    for i=1:length(trialData)
        trialData(i).dayLabel = d;
        trialData(i).newTrial = trialData(i).trial + trialNumBegin;
    end

    trialDataAll = cat(1, trialDataAll, trialData); clear trialData
    behavDataAll = cat(2, behavDataAll, behavData); clear behavData
    trialNumBegin = trialNumBegin + max(trialNum);
end; clear taskfile date dataFolder
trialData = trialDataAll; clear trialDataAll
behavData = behavDataAll; clear behavDataAll

%% get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
trialStatusLabels = [behavData.trialStatusLabel];
delayTimes = [behavData.delayTime]; delayTimes = round(delayTimes / 50) * 50;
reactionTimes = [behavData.reactionTime];
peakSpeeds = [behavData.peakSpeed];
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
DiffStyle = ["-", ":"];
DelayTimes = 300:50:1050; nDelayTimes = length(DelayTimes);

%% SUCCESS RATE
failedLabels = trialStatusLabels ~= 1;
VS_reward_size_success_rate(failedLabels, rewardLabels, difficultyLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
VS_delay_time_success_rate(failedLabels, rewardLabels, difficultyLabels, delayTimes, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
VS_direction_success_rate(failedLabels, rewardLabels, difficultyLabels, directionLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
VS_trial(trialNums, 100 * int16(~failedLabels), rewardLabels, difficultyLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in preparation rate
delayFailure = trialStatusLabels == -11 | trialStatusLabels == -12 | trialStatusLabels == 0;
VS_reward_size_success_rate(delayFailure, rewardLabels, difficultyLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderFailureMode+"SRInPreparation", "Ylim", [60 100])
VS_delay_time_success_rate(delayFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderFailureMode+"SRInPreparation", "Ylim", [60 100])
VS_direction_success_rate(delayFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderFailureMode+"SRInPreparation", "Ylim", [60 100])
VS_trial(trialNums, 100 * int16(~delayFailure), rewardLabels, difficultyLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderFailureMode+"SRInPreparation", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in reaching rate
reachFailure = ismember(trialStatusLabels, [-22, -23, -24]);
VS_reward_size_success_rate(reachFailure, rewardLabels, difficultyLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderFailureMode+"SRInReaching", "Ylim", [60 100])
VS_delay_time_success_rate(reachFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderFailureMode+"SRInReaching", "Ylim", [60 100])
VS_direction_success_rate(reachFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderFailureMode+"SRInReaching", "Ylim", [60 100])
VS_trial(trialNums, 100 * int16(~reachFailure), rewardLabels, difficultyLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderFailureMode+"SRInReaching", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in holding rate
holdFailure = ismember(trialStatusLabels, [-31, -32, -33, -34, -35]);
VS_reward_size_success_rate(holdFailure, rewardLabels, difficultyLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderFailureMode+"SRInHolding", "Ylim", [60 100])
VS_delay_time_success_rate(holdFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in holding rate (%)", "OutputFolder", outputFolderFailureMode+"SRInHolding", "Ylim", [60 100])
VS_direction_success_rate(holdFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderFailureMode+"SRInHolding", "Ylim", [60 100])
VS_trial(trialNums, 100 * int16(~holdFailure), rewardLabels, difficultyLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderFailureMode+"SRInHolding", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% REACTION TIME
VS_reward_size(reactionTimes, rewardLabels, difficultyLabels, "Reaction Time (ms)", outputFolderVigor+"ReactionTime")
VS_delay_time(reactionTimes, rewardLabels, difficultyLabels, delayTimes, "Reaction Time (ms)", outputFolderVigor+"ReactionTime")
VS_direction(reactionTimes, rewardLabels, difficultyLabels, directionLabels, "Reaction Time (ms)", outputFolderVigor+"ReactionTime")
VS_trial(trialNums, reactionTimes, rewardLabels, difficultyLabels, "Label", "Reaction Time (ms)", "OutputFolder", outputFolderVigor+"ReactionTime", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% PEAK SPEED
VS_reward_size(peakSpeeds, rewardLabels, difficultyLabels, "Peak Speed (m/s)", outputFolderVigor+"PeakSpeed")
VS_delay_time(peakSpeeds, rewardLabels, difficultyLabels, delayTimes, "Peak Speed (m/s)", outputFolderVigor+"PeakSpeed")
VS_direction(peakSpeeds, rewardLabels, difficultyLabels, directionLabels, "Peak Speed (m/s)", outputFolderVigor+"PeakSpeed")
VS_trial(trialNums, peakSpeeds, rewardLabels, difficultyLabels, "Label", "Peak Speed (m/s)", "OutputFolder", outputFolderVigor+"PeakSpeed", ...
        "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)