close all; clear all;
addpath(genpath("./utils/plotting/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
ndates = length(dates);
axisOutputFileName = "../interim/GC_rewardAxis_allDay.mat";
outputFolder = "../results/202306w3-results-summary-all/neuralData_GoCue/rewardAxis_vs_parameters/"; makeDir(outputFolder);
outputFolder1DProjection = "../results/202306w3-results-summary-all/neuralData_GoCue/rewardAxisProjection/"; makeDir(outputFolder1DProjection);
analysisBin = (50:250); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(400:600)

trialDataAll = struct.empty;
trialNumBegin = 0;
for d = 1:ndates
    date = dates(d);
    load("../data/processed/stitched-whole/2022"+date+"_GC_200_600_3456.mat");
%     load("../data/processed/stitched-whole/2022"+date+"_TO_200_700_3.mat");
%     load("../data/processed/stitched-whole/2022"+date+"_HT_200_375_3456.mat");
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
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
trialStatusLabels = [trialData.trialStatus];
delayTimes = [trialData.delayTime]; delayTimes = round(delayTimes / 50) * 50;

% Get data and avg within dir x rew
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);
neuralData_byParameters_mean = nan(ndirections,nrewards,ndifficulties,nneurons);
for l = 1:ndirections
    for r = 1:nrewards
        for j = 1:ndifficulties
            curInds = directionLabels==directions(l) & rewardLabels==rewards(r) & difficultyLabels==difficulties(j);
            neuralData_byParameters_mean(l,r,j,:) = mean(neuralData(curInds,:));
        end; clear j
    end; clear r
end; clear l

% reward axis of both
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 3]));
% meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zR(3, 1) < zR(1, 1)
    wR(:, 1) = -wR(:, 1);
end

% reward axis of Tiny
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 1, :));
[wRTiny,zRTiny,eigVls_RTiny] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRTiny(3, 1) < zRTiny(1, 1)
    wRTiny(:, 1) = -wRTiny(:, 1);
end

% reward axis of Huge
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wRHuge,zRHuge,eigVls_RHuge] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRHuge(3, 1) < zRHuge(1, 1)
    wRHuge(:, 1) = -wRHuge(:, 1);
end

% difficulty axis
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 2]));
[wDif,zDif,eigVls_Dif]  = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndifficulties-1);
if zDif(2, 1) < zDif(1, 1)
    wDif(:, 1) = -wDif(:, 1);
end

% R-Diff axis
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = reshape(meanNeuralData, [nrewards*ndifficulties, nneurons]);
[wRDiff,zRDiff,eigVls_RDiff] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards*ndifficulties-1);
if zRDiff(3, 1) < zRDiff(1, 1)
    wRDiff(:, 1) = -wRDiff(:, 1);
end
if zRDiff(4, 2) < zRDiff(1, 2)
    wRDiff(:, 2) = -wRDiff(:, 2);
end

save(axisOutputFileName, "wR", "wRTiny", "wRHuge", "wDif", "wRDiff", "eigVls_R", "eigVls_RTiny", "eigVls_RHuge", "eigVls_Dif", "eigVls_RDiff");

%% plot reward PC1 and difficulty PC1 by 2D scatter plot
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRDiff(:,1:2);
figure;
hold on;
for i=1:ndifficulties
    for j=1:nrewards
        curInds = rewardLabels==rewards(j) & difficultyLabels==difficulties(i);
        X(i,j) = mean(neuralData_onDPC(curInds, 1));
        Xerr(i,j) = std(neuralData_onDPC(curInds, 1))/sqrt(sum(curInds));
        Y(i,j) = mean(neuralData_onDPC(curInds, 2));
        Yerr(i,j) = std(neuralData_onDPC(curInds, 2))/sqrt(sum(curInds));
    end
    h(i) = errorbar(X(i,:), Y(i,:), Yerr(i,:), Yerr(i,:), Xerr(i,:), Xerr(i,:), "Color", diffColors(i,:), 'LineWidth', 2.5); hold on;
    scatter(X(i,:), Y(i,:), 100, rewColors, "filled"); hold on;
end
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
ylim([-30 20]); xlim([-30 20]);
xlabel("RDiff Axis: PC1 (Hz): "+round(eigVls_RDiff(1)/sum(eigVls_RDiff)*100, 2)+"%");
ylabel("RDiff Axis: PC2 (Hz): "+round(eigVls_RDiff(2)/sum(eigVls_RDiff)*100, 2)+"%");
legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
saveas(gcf, outputFolder+"2Dprojection-all.jpg");

%% plot reward PC1 and difficulty PC1
rAxes = ["Tiny", "Huge", "All"];
for k=1:3 % 1:2
    Y = zeros(3,2); % reward x difficulty
    Yerr = zeros(3,2);
    if k==1
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRTiny(:,1);
    elseif k==2
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRHuge(:,1);
    else
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
    end

    for i=1:nrewards
        for j=1:ndifficulties
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j);
            Y(i,j) = mean(neuralData_onDPC(curInds));
            Yerr(i,j) = std(neuralData_onDPC(curInds))/sqrt(sum(curInds));
        end
    end
    figure;
    for j=1:ndifficulties
        errorbar(1:nrewards, Y(:,j), Yerr(:,j), "Color", diffColors(j,:), 'LineWidth', 2); hold on;
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xticks([1 2 3]); xticklabels(["S", "M", "L"]); xlim([0.7 3.3]); legend(["Tiny", "Huge"], Location="best")
    set(gcf,'position',[100,100,400,650]);
    ylabel("Reward Axis (Hz)");
    saveas(gcf, outputFolder1DProjection+"Axis-in-"+rAxes(k)+".jpg"); close all;
end


neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);


%% Reward Axis
VS_reward_size(neuralData_onDPC', rewardLabels, difficultyLabels, "Reward Axis (spikes/s)", outputFolder+"RewardAxis")
VS_delay_time(neuralData_onDPC', rewardLabels, difficultyLabels, delayTimes, "Reward Axis (spikes/s)", outputFolder+"RewardAxis")
VS_direction(neuralData_onDPC', rewardLabels, difficultyLabels, directionLabels, "Reward Axis (spikes/s)", outputFolder+"RewardAxis")
VS_trial(trialNums, neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolder+"RewardAxis", ...
            "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1650)