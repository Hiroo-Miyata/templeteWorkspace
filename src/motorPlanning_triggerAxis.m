close all;clear all;
addpath(genpath("./utils/plotting/"));
addpath(genpath("./utils/processing/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
Yall = zeros(length(dates), 3, 2, 8);
Yall2 = zeros(length(dates), 3, 2, 8);
for d = 1:length(dates)
date = dates(d);
outputFolder = "../results/202306w3-results-summary-all/neuralData_GoCue/motorPlanning_triggerAxis/"; makeDir(outputFolder);
% load trialData
load("../data/processed/non-stitched/2022"+date+"_GC_200_600_3456.mat");
taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
load(taskfile.folder +"/"+ taskfile.name);
trialNum = [trialData.trial];

behavData = behavData(ismember([behavData.trial], trialNum));
for i = 1:length(trialData)
    trialData(i).reactionTime = behavData(i).reactionTime;
    trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
end


% Get labels
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
reactionTimes = [trialData.reactionTime];

curInds = ~isnan(reactionTimes) & reactionTimes > 200 & reactionTimes < 450;
trialData = trialData(curInds);


directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
reactionTimes = [trialData.reactionTime];
movementOnset = round(reactionTimes) + 200; % binWindow starts 200 ms before GO Cue

% get neural data
neuralActivity = cat(3, trialData(:).firingRates);
[neuralActivity, nneurons, ntimebins, ntrials] = kernelSmooth(neuralActivity);

% Requirements of SVM
% The SVM should be optimized to distinguish between points occurring 360-180 ms prior to movement onset and points occurring between 120 ms prior to movement onset and 60 ms after movement onset.
% Performance should be assessed on a test set of 10% of trials using leave-one-out cross-validation.
% The SVM should identify the dimension along which neural activity changes prior to a reach.
% prepare dataset


% separate -1 and 1
beforeTimes = -360:10:-240;
nbeforeT = length(beforeTimes);
afterTimes = -60:10:60;
nafterT = length(afterTimes);
beforeReachData = zeros(nneurons, nbeforeT, ntrials);
afterReachData = zeros(nneurons, nafterT, ntrials);
for i = 1:ntrials
    beforeReachData(:, :, i) = neuralActivity(:, movementOnset(i)+beforeTimes, i);
    afterReachData(:, :, i) = neuralActivity(:, movementOnset(i)+afterTimes, i);
end
beforeReachData = reshape(beforeReachData, nneurons, nbeforeT*ntrials);
afterReachData = reshape(afterReachData, nneurons, nafterT*ntrials);
beforeReachData = beforeReachData.';
afterReachData = afterReachData.';
% normalize in each neuron across all trials
Zmean = mean([beforeReachData; afterReachData]);
Zstd = std([beforeReachData; afterReachData]);
beforeReachData = (beforeReachData - Zmean) ./ Zstd;
afterReachData = (afterReachData - Zmean) ./ Zstd;
% concatenate
data = [beforeReachData; afterReachData];
% label
label = [ones(size(beforeReachData, 1), 1); -ones(size(afterReachData, 1), 1)];

% save("../data/processed/SVM/SVMData_"+date+".mat", "data", "label");
%% LDA
% requirement: Performance should be assessed on a test set of 10% of trials using leave-one-out cross-validation.
% split data into training and test set
% shuffle the dataset before the cross validation and keep the distribution of labels same across folds

kfolds = 10;
fold_size = floor(length(label)/kfolds);

% shuffle
rng(1);
shuffledIndex = randperm(length(label));
shuffledData = data(shuffledIndex, :);
shuffledLabel = label(shuffledIndex);

% cross validation
% store trigger dimension which is orthogonal to the hyperplane
% accuracy = zeros(kfolds, 1);
% hyperPlane = zeros(nneurons, kfolds);
% for i = 1:kfolds
%     % split data
%     testIndex = (i-1)*fold_size+1:i*fold_size;
%     trainIndex = setdiff(1:length(label), testIndex);
%     trainData = shuffledData(trainIndex, :);
%     testData = shuffledData(testIndex, :);
%     trainLabel = shuffledLabel(trainIndex);
%     testLabel = shuffledLabel(testIndex);
    
%     % train LDA
%     Mdl = fitcdiscr(trainData, trainLabel);
%     % test LDA
%     predictedLabel = predict(Mdl, testData);
%     accuracy(i) = sum(predictedLabel == testLabel) / length(testLabel);
% end

% disp(mean(accuracy));


%% LDA in all trials

% train LDA
Mdl = fitcdiscr(data, label);
% plot the data projected on the vector calulated by LDA
% get the vector
wTrigger = Mdl.Coeffs(1,2).Linear;
save("../interim/GC_TriggerAxis_"+date+".mat", "wTrigger");
% project the data
projectedDataBefore = beforeReachData * wTrigger;
projectedDataAfter = afterReachData * wTrigger;

% reshape the projectedData
beforeOnLDA = reshape(projectedDataBefore, nbeforeT, ntrials);
afterOnLDA = reshape(projectedDataAfter, nafterT, ntrials);

% plot
% before color is red, after color is blue
% figure;
% hold on;
% for i = 1:ntrials 
%     % scatter plot
%     scatter(beforeTimes, beforeOnLDA(:, i), 'filled', 'MarkerFaceColor', 'r');
%     scatter(afterTimes, afterOnLDA(:, i), 'filled', 'MarkerFaceColor', 'b');
% end
% title("Data projected on the trigger dimension");
% xlabel("Time (ms)");
% ylabel("Projection on the trigger dimension");
% saveas(gcf, outputFolder+"LDA_projection_"+date+".png");

%% plot the trajectory of the trigger dimension in each difficulty
[nneurons, binWindow, ntrials] = size(neuralActivity);

figure; hold on;
for i=1:ndifficulties
    meanFR = mean(neuralActivity(:, :, difficultyLabels==difficulties(i)), 3);
    meanFR = squeeze(meanFR);
    meanFR = meanFR.';
    % project the data
    projectedData = meanFR * wTrigger;
    % plot
    plot(-200:600, projectedData, 'LineWidth', 3, "Color", diffColors(i, :));
end
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
xlabel("Time (ms)");
xticks(-200:200:600);
xticklabels(["-200", "GC", "200", "400", "600"]);
ylabel("Projection on the triggerAxis (a.u.)");
legend("Tiny", "Huge"); set(gcf,'position',[0,0,550,550]); 
saveas(gcf, outputFolder+"LDA_projection_difficulty_"+date+".png");

%% calculate the correlation between the value on trigger dimension at GoCue and the reaction time
analysisBin = (201:201); % Pay Attention!! GC=(50:250), TO=(150:350)
rTcorr = zeros(nrewards, ndifficulties, ndirections);
for i = 1:nrewards
    for j = 1:ndifficulties
        for k = 1:ndirections
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
            meanFR = mean(neuralActivity(:, analysisBin, :), 2);
            meanFR = squeeze(meanFR);
            meanFR = meanFR.';
            TriggerScore = meanFR * wTrigger;
            rTcorr(i, j, k) = corr(TriggerScore(curInds), reactionTimes(curInds).');
        end
    end
end
Yall(d, :, :, :) = rTcorr; 


%% get the correlation between the time over threshold and reaction time
rTcorr = zeros(nrewards, ndifficulties, ndirections);
for i = 1:nrewards
    for j = 1:ndifficulties
        for k = 1:ndirections
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
            curIdxs = find(curInds);
            triggerTime = zeros(length(curIdxs), 1);
            for l = 1:length(curIdxs)
                neuralData = squeeze(neuralActivity(:, :, curIdxs(l)));
                neuralData = neuralData.';
                score = neuralData * wTrigger;
                % find the time over threshold
                threshold = (max(score) + min(score)) / 2;
                % find the time over threshold
                triggerTime(l) = find(score > threshold, 1);
            end
            rTcorr(i, j, k) = corr(triggerTime, reactionTimes(curInds).');
        end
    end
end
Yall2(d, :, :, :) = rTcorr; 
close all;
end

%% plotting
corrHistogram(Yall(:), "Label", "Correlation Coefficient", "Title", "RT vs TriggerAxis-at-GC", "OutputFolder", outputFolder+"all-at-GC")
corrHistogram(Yall2(:), "Label", "Correlation Coefficient", "Title", "RT vs TriggerAxis-at-Threshold", "OutputFolder", outputFolder+"all-at-Threshold")