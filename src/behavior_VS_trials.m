close all; clear all;
rootDir = "../";
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
for d = 1:length(dates)
date = dates(d);
outputFolder = "../results/202302w1/RTTrajectory/";
%% data loading
dataFolder = rootDir + "data/preprocessed/2022"+date+"/";
taskfile = dir(dataFolder+ "*" +date+ "*_TSK_*.mat");
load(taskfile.folder +"/"+ taskfile.name);
kinfile = dir(dataFolder+ "*" +date+ "*_KIN_*.mat");
load(kinfile.folder +"/"+ kinfile.name);
trialData = arrayfun(@(x) setfield(trialData(x), "handKinematics", kinematicData(x)),1:length(trialData));

validTrialData = struct.empty(0);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};

for i=(1:length(trialData))
    stateTransition = trialData(i).stateTable;
    if all(ismember([3 4], stateTransition(1,:))) == 1
        validTrialData = cat(1, validTrialData, trialData(i));
    elseif all(ismember([3 11], stateTransition(1,:))) == 1
        validTrialData = cat(1, validTrialData, trialData(i));
    else
        
    end
end
tasktmp = taskInfo.taskParams;
tasktmp.subjectName = taskInfo.subjectName;
updatedTrialData = p5_updateBehavior(validTrialData, tasktmp);
clear trialData;

directionLabels = [updatedTrialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [updatedTrialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, updatedTrialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

% get reaction time
rT = [updatedTrialData.t4_reactionTime];
GoCueTime = [updatedTrialData.t3_goCueTime];
rT = rT - GoCueTime;
pS = [updatedTrialData.peakSpeed];
curInds = ~isnan(rT) & rT > 200 & rT < 450;
rT = rT(curInds);
dLs = directionLabels(curInds);
diffLs = difficultyLabels(curInds);
rewLs = rewardLabels(curInds);

%% plot the smoothed trajectory of RT in each reward and difficulty condition
figure;
for i = 1:nrewards
    for j = 1:ndifficulties
        curInds = rewLs == rewards(i) & diffLs == difficulties(j);
        curRT = rT(curInds);
        smoothedRT = movmean(curRT, 50);
        h(j) = plot(smoothedRT, 'Color', diffColors(j, :), 'LineWidth', 2); hold on;
        scatter([1 length(smoothedRT)], [smoothedRT(1) smoothedRT(end)], 100, 'filled', 'MarkerFaceColor', rewColors(i, :)); hold on;
    end
end
legend(h, {"Tiny", "Huge"});
title("Reaction Time over time");
xlabel("Trial Number");
ylabel("Reaction Time (ms)");
saveas(gcf, outputFolder + "RT_" + date + ".jpg");

%% Peak Speed
pS = [updatedTrialData.peakSpeed];
curInds = ~isnan(pS);
pS = pS(curInds);
dLs = directionLabels(curInds);
diffLs = difficultyLabels(curInds);
rewLs = rewardLabels(curInds);

%% plot the smoothed trajectory of RT in each reward and difficulty condition
VS_trial(trialNums, Y, rewardLabels, difficultyLabels, options)
figure;
for i = 1:nrewards
    for j = 1:ndifficulties
        curInds = rewLs == rewards(i) & diffLs == difficulties(j);
        curPS = pS(curInds);
        smoothedPS = movmean(curPS, 50);
        h(j) = plot(smoothedPS, 'Color', diffColors(j, :), 'LineWidth', 2); hold on;
        scatter([1 length(smoothedPS)], [smoothedPS(1) smoothedPS(end)], 100, 'filled', 'MarkerFaceColor', rewColors(i, :)); hold on;
    end
end
legend(h, {"Tiny", "Huge"});
title("Peak Speed over time");
xlabel("Trial Number");
ylabel("Peak Speed (m/s)");
saveas(gcf, outputFolder + "PS_" + date + ".jpg");

%% delay success rate vs time


close all;
end
