close all;clear all;
addpath(genpath("./utils/plotting/"));
addpath(genpath("./utils/processing/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
Yall = zeros(length(dates), 3, 2, 8);
Yall2 = zeros(length(dates), 3, 2, 8);
for d = 1:length(dates)
date = dates(d);
outputFolder = "../results/202306w3-results-summary-all/neuralData_GoCue/motorPlanning_DynamicPCs/";
% load trialData
load("../data/processed/non-stitched/2022"+date+"_GC_200_600_3456.mat");
taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
load(taskfile.folder +"/"+ taskfile.name);
trialNum = [trialData.trial];

behavData = behavData(ismember([behavData.trial], trialNum));
for i = 1:length(trialData)
    trialData(i).reactionTime = behavData(i).zReactionTime;
    trialData(i).peakSpeed = behavData(i).zPeakSpeed;
end

% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
reactionTimes = [trialData.reactionTime];
peakSpeeds = [trialData.peakSpeed];

% Get mean PSTH in each condition
neuralActivity = cat(3, trialData(:).firingRates);
[neuralActivity, nneurons, ntimebins, ntrials] = kernelSmooth(neuralActivity);

% 1. get the average firing rate of more than 700 ms delay period data
% 2. get the average firing rate of delay successed trials
% 3. get the average firing rate of both delay successed and delay failed trials

timeBin = 501;

% 0. get the average firing rate of all trials
meanPSTHafterGC = zeros(ndirections, nrewards, ndifficulties, timeBin, nneurons);
for i = 1:ndirections
    for j = 1:nrewards
        for k = 1:ndifficulties
            curInds = directionLabels==directions(i) & rewardLabels==rewards(j) & difficultyLabels==difficulties(k);
            PSTHeach = neuralActivity(:, 1:timeBin, curInds);
            meanPSTHafterGC(i, j, k, :, :) = squeeze(mean(PSTHeach, 3))';
        end
    end
end


meanPSTHafterGC = mean(meanPSTHafterGC, [2 3]);
meanPSTHafterGC = reshape(meanPSTHafterGC, [ndirections*timeBin, nneurons]); %[nrewards*ndifficulties*timeBin, nneurons]

meanPSTHafterGC = meanPSTHafterGC-mean(meanPSTHafterGC, 1);
[wDynamic,zDynamic,eigVls_Dynamic] = pca(meanPSTHafterGC);

save("../interim/GC_dynamicsAxis_"+date+".mat", "wDynamic", "zDynamic", "eigVls_Dynamic")


trajectoryOnPCs = nan(ntrials, timeBin, 3);
for i=1:ntrials
    PSTH = squeeze(neuralActivity(:, 1:timeBin, i)); % nneurons * ntimeBin
    trajectoryOnPCs(i, :, :) = PSTH' * wDynamic(:, 1:3);
end; clear i PSTH

eachTrajectory = zeros(ndirections, nrewards, ndifficulties, timeBin, 3);
for i = 1:ndirections
    for j = 1:nrewards
        for k = 1:ndifficulties
            curInds = directionLabels==directions(i) & rewardLabels==rewards(j) & difficultyLabels==difficulties(k);
            eachTrajectory(i, j, k, :, :) = mean(trajectoryOnPCs(curInds, :, :), 1);
            % eachTrajectory_std(i, j, k, :, :) = std(trajectoryOnPCs, 1) / sqrt(sum(~isnan(trajectoryOnPCs(curInds, 1, 1))));
        end
    end
end

eachTrajectory = squeeze(mean(eachTrajectory, 2));
% eachTrajectory_std = squeeze(mean(eachTrajectory_std, 2));
% show the animation of the trajectory of the first two principal components
pc = 2;
figure; hold on;
for i = 1:ndirections
    for k = 1:ndifficulties
        % show eachTrajectory by animation
        timeBin = 1:401;
        X = squeeze(eachTrajectory(i, k, timeBin, 1));
        Y = squeeze(eachTrajectory(i, k, timeBin, pc));
        % Xstd = squeeze(eachTrajectory_std(i, k, timeBin, 1));
        % Ystd = squeeze(eachTrajectory_std(i, k, timeBin, pc));
        % plot the line
        h(k) = plot(X, Y, 'Color', diffColors(k, :), 'LineWidth', 2);
        % patch([X+Xstd; flip(X-Xstd)], [Y+Ystd; flip(Y-Ystd)], rewColors(j, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none')        % plot the error range
        plot(X(200), Y(200), 'o', 'Color', diffColors(k, :), 'MarkerSize', 12, 'MarkerFaceColor', diffColors(k, :))
    end
end
% xlabel and ylabel add the percentage of the variance explained by the first two principal components
xlabel("PC1 (Hz) "+round(eigVls_Dynamic(1)*100/sum(eigVls_Dynamic))+"%");
ylabel("PC"+pc+" (Hz) "+round(eigVls_Dynamic(pc)*100/sum(eigVls_Dynamic))+"%");
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold'); legend(h, ["Tiny", "Huge"], Location="best");
set(gcf,'position',[0,0,550,550]);
saveas(gcf, outputFolder+"Around-GC-PC"+pc+"-"+date+".jpg");
close all
fprintf("PC1: %d, PC2: %d, PC3: %d, PC4: %d \n", round(eigVls_Dynamic(1:4)*100/sum(eigVls_Dynamic)))


rTcorr = zeros(nrewards, ndifficulties, ndirections);
meanDynamicPC = zeros(nrewards, ndifficulties, ndirections);
TriggerScore = squeeze(neuralActivity(:, 201, :))' * wDynamic(:, 1);
RTInds = ~isnan(reactionTimes);

for i = 1:nrewards
    for k = 1:ndirections
        curInds = rewardLabels==rewards(i) & directionLabels==directions(k);
        zTriggerScore = nan(size(TriggerScore));
        zTriggerScore(curInds) = zscore(TriggerScore(curInds));
        for j = 1:ndifficulties
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
            curInds = curInds & RTInds;
            rTcorr(i, j, k) = corr(TriggerScore(curInds), reactionTimes(curInds).');
            meanDynamicPC(i,j,k) = mean(TriggerScore(curInds));
        end
    end
end
Yall(d, :, :, :) = rTcorr; 
Yall2(d, :, :, :) = meanDynamicPC;
% corrHistogram(rTcorr(:), "Label", "Correlation Coefficient", "Title", "RT vs Dynamic PC1", "OutputFolder", outputFolder+date)
end
% corrHistogram(Yall(:), "Label", "Correlation Coefficient", "Title", "RT vs Dynamic PC1", "OutputFolder", outputFolder+"all")

x = squeeze(Yall2(:, :, 1, :));
y = squeeze(Yall2(:, :, 2, :));
rankTest2d(x(:), y(:), "XLabel", "Tiny Position", "YLabel", "Huge Position", "OutputFolder", outputFolder+"position-at-GC")
corrHistogram(y(:) - x(:), "Label", "Huge Position - Tiny Position", "Title", "Position shift at GC", "OutputFolder", outputFolder+"position-at-GC", ...
                "Xlim", "none")