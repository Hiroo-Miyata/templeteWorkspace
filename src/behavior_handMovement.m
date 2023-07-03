close all; clear all;addpath('./utils/plotting/');
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
fileName = "all";
outputFolder = "../results/202306w3-results-summary-"+ fileName + "/behavior_hand-movement/";
makeDir(outputFolder)

%% conbine across days
trialDataAll = struct.empty;
trialNumBegin = 0;
for d = 1:length(dates)
    date = dates(d);
    load("../data/processed/stitched/2022"+date+"_GC_200_700_34567.mat");
    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    trialNum = [trialData.trial];

    behavData = behavData(ismember([behavData.trial], trialNum));
    for i = 1:length(trialData)
        trialData(i).reactionTime = behavData(i).t4_reactionTime - behavData(i).t3_goCueTime;
        trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
        trialData(i).dayLabel = d;
        trialData(i).newTrial = trialData(i).trial + trialNumBegin;
    end
    trialDataAll = cat(1, trialDataAll, trialData);
    trialNumBegin = trialNumBegin + max(trialNum);
end
trialData = trialDataAll;


%% get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
reactionTimes = [trialData.reactionTime];
delayTimes = [trialData.delayTime];
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
DiffStyle = ["-", ":"];
DelayTimes = 300:50:1050; nDelayTimes = length(DelayTimes);


% calculate the error of hand movement in each trial
% the error rate until the target acquisition
movementError = nan(1, length(trialData));
for i = 1:length(trialData)
    handPositions = trialData(i).handKinematics.rotatedPosition;
    movementOnset = reactionTimes(i) + 200;

    if ~isnan(movementOnset)
        stateTransition = trialData(i).stateTable;
        TargetAcquisition = stateTransition(2, find(stateTransition(1, :)==6)) - stateTransition(2, find(stateTransition(1, :)==4)) + 200;
        if TargetAcquisition > 900
            disp("Error occur");
        end
        movement= handPositions(int16(movementOnset):TargetAcquisition, :);
        movementError(i) = 0;
        for j = 2:size(movement, 1)
            deltaX = abs(movement(j, 1) - movement(j-1, 1));
            if (movement(j, 2) < 0 & movement(j-1, 2) < 0) & (movement(j, 2) > 0 & movement(j-1, 2) > 0)
                deltaY = (abs(movement(j, 2)) + abs(movement(j-1, 2))) * 0.25;
            else
                deltaY = abs(movement(j, 2) + movement(j-1, 2)) * 0.5;
            end
            deltaMove = deltaX * deltaY;
            movementError(i) = movementError(i) + deltaMove;
        end
        if movementError(i) > 500
            movementError(i) = nan;
        end
    end
end

%% Movement Error
VS_reward_size(movementError, rewardLabels, difficultyLabels, "Movement Error (a.u.)", outputFolder+"MovementError")
VS_delay_time(movementError, rewardLabels, difficultyLabels, delayTimes, "Movement Error (a.u.)", outputFolder+"MovementError")
VS_direction(movementError, rewardLabels, difficultyLabels, directionLabels, "Movement Error (a.u.)", outputFolder+"MovementError")
VS_trial(trialNums, movementError, rewardLabels, difficultyLabels, "Label", "Movement Error (a.u.)", "OutputFolder", outputFolder+"MovementError", ...
            "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1700)