close all;clear all;
addpath(genpath("./utils/plotting/"));
addpath(genpath("./utils/processing/"));

% main function 
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
outputFolder = "../results/202306w2/vigorAxisAnalaysis_RT/Ridge";
outputFolderSub = "../results/202306w2/vigorAxis_RT_regressionResult_Ridge/";
outputFolderRT = "../results/202306w2/RT_PDF_eachDay/";
outputFolderAngle = "../results/202306w2/vigorAxis-Angle/Ridge";

trialDataAll = struct.empty;
for d = 1:length(dates)
    date = dates(d);
    load("../data/processed/stitched/2022"+date+"_GC_200_600_3456.mat");
    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    trialNum = [trialData.trial];
    trialData = trialData(trialNum > 24);
    trialNum = [trialData.trial];
    behavData = behavData(ismember([behavData.trial], trialNum));
    for i = 1:length(trialData)
        % trialData(i).reactionTime = behavData(i).t4_reactionTime - behavData(i).t3_goCueTime;
        % trialData(i).peakSpeed = behavData(i).peakSpeed;
        trialData(i).reactionTime = behavData(i).zReactionTime;
        trialData(i).peakSpeed = behavData(i).zPeakSpeed;
        trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
        trialData(i).trialNew = trialData(i).trial + d*10000;
    end
    reactionTimes = [trialData.reactionTime];
    trialData = trialData(~isnan(reactionTimes));

    trialDataAll = cat(1, trialDataAll, trialData);

    % in each trialData set, see the reaction time projection. The histgram show the PDF of the reaction time
    reactionTimes = [trialData.reactionTime];
    [counts, edges] = histcounts(reactionTimes, 50, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end))/2;
    figure;
    bar(centers, counts, 'FaceColor', [0.7 0.7 0.7]);  % change color as needed
    xlabel('Reaction Time (ms)');
    ylabel('PDF');
    title('Histogram with PDF on Y-axis');
    saveas(gcf, outputFolderRT + "RT_PDF_zscored_" + date + ".png");


end
trialData = trialDataAll; clear trialDataAll

% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trialNew];
reactionTimes = [trialData.reactionTime];
peakSpeeds = [trialData.peakSpeed];

% Get data and avg within dir x rew
analysisBin = (50:250); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(100:300)
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);

% regression sub group
% 各reward/difficulty <= どのreward/difficultyにおいても、リアクションタイムが分散を大きくする
% 全体 <= 各方向でoptimal subspaceが存在し、それがreaction timeと関係している

%% find vigor axis for each direction
% by using findVigorAxisCV function
% [W, W0, r, best_lambda, FitInfo] = findVigorAxisCV(neuralData, vigorData, kfolds)

vigorAxes_eachDir = zeros(ndirections, nneurons);
vigorIntercepts_eachDir = zeros(ndirections, 1);
vigorAxes_eachDir_Diff_rew = zeros(ndirections, ndifficulties, nrewards, nneurons);
vigorIntercepts_eachDir_Diff_rew = zeros(ndirections, ndifficulties, nrewards);

for i = 1:ndirections
    % find vigor axis for each direction
    directionIdx = directionLabels == directions(i);
    vigorData = reactionTimes(directionIdx);
    [W, W0, r, best_lambda, FitInfo] = findVigorAxisCV(neuralData(directionIdx, :), vigorData', 5, outputFolderSub + "vigorAxis_eachDir_" + directions(i));
    vigorAxes_eachDir(i, :) = W;
    vigorIntercepts_eachDir(i) = W0;
    
    % find vigor axis for each direction, difficulty, reward
    for j = 1:ndifficulties
        for k = 1:nrewards
            idx = directionIdx & difficultyLabels == difficulties(j) & rewardLabels == rewards(k);
            vigorData = reactionTimes(idx);
            [W, W0, r, best_lambda, FitInfo] = findVigorAxisCV(neuralData(idx, :), vigorData', 5, ...
                                            outputFolderSub + "vigorAxis_eachDir_Diff_rew_" + directions(i) + "_" + difficulties(j) + "_" + rewards(k));
            vigorAxes_eachDir_Diff_rew(i, j, k, :) = W;
            vigorIntercepts_eachDir_Diff_rew(i, j, k) = W0;
        end
    end
end

%% plot vigor axis onto reward-difficulty space for each direction
neuralData_byParameters_mean = nan(ndirections,nrewards,ndifficulties,nneurons);
for l = 1:ndirections
    for r = 1:nrewards
        for j = 1:ndifficulties
            curInds = directionLabels==directions(l) & rewardLabels==rewards(r) & difficultyLabels==difficulties(j);
            neuralData_byParameters_mean(l,r,j,:) = mean(neuralData(curInds,:));
        end; clear j
    end; clear r
end; clear l

RDiffAxes = zeros(ndirections, nneurons, 2);
angles = zeros(ndirections, 2);
for i = 1:ndirections
    meanNeuralData = squeeze(neuralData_byParameters_mean(i, :, :, :));
    % meanNeuralData = squeeze(mean(neuralData_byParameters_mean(:, :, :, :), 1));
    meanNeuralData = reshape(meanNeuralData, [nrewards*ndifficulties, nneurons]);
    [wRDiff,zRDiff,eigVls_RDiff] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards*ndifficulties-1);
    if zRDiff(3, 1) < zRDiff(1, 1)
        wRDiff(:, 1) = -wRDiff(:, 1);
    end
    if zRDiff(4, 2) < zRDiff(1, 2)
        wRDiff(:, 2) = -wRDiff(:, 2);
    end
    RDiffAxes(i, :, :) = wRDiff(:, 1:2);

    meanNeuralData = squeeze(neuralData_byParameters_mean(i, :, :, :));
    meanNeuralData = squeeze(mean(meanNeuralData, 1));
    meanNeuralData = reshape(meanNeuralData, [ndifficulties, nneurons]);
    [wDif,zDif,eigVls_Dif]  = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndifficulties-1);
    if zDif(2, 1) < zDif(1, 1)
        wDif(:, 1) = -wDif(:, 1);
    end

    % plot vigor axis
    % X axis: wRDiff(:, 1), Y axis: wRDiff(:, 2)
    % plot the mean projection in each reward and difficulty condition
    % plot the vigor axis of each reward and difficulty condition
    % which is arrow from mean to mean + 10
    % add text next to the arrow to indicate the angle between the vigor axis
    % and the wRDiff(:, 1:2)
    % Also project the vigor axis of each direction onto the RDiffAxes
    % the arrow should be from the xlim and ylim of the RDiffAxes
    % the text should be next to the arrow
    % the text should indicate the angle between the vigor axis and the
    % RDiffAxes
    figure('Name', "vigor axis for direction " + directions(i))
    hold on
    for j = 1:ndifficulties
        for k = 1:nrewards
            curInds = directionLabels==directions(i) & rewardLabels==rewards(k) & difficultyLabels==difficulties(j);
            curNeuralData = mean(neuralData(curInds, :), 1);
            curNeuralData_onPC = curNeuralData * wRDiff(:, 1:2);
            curVigorAxis = vigorAxes_eachDir_Diff_rew(i, j, k, :);
            curVigorAxis = squeeze(curVigorAxis);
            projVigorAxis(1) = dot(curVigorAxis, squeeze(RDiffAxes(i, :, 1))) / norm(curVigorAxis) / norm(squeeze(RDiffAxes(i, :, 1)));
            projVigorAxis(2) = dot(curVigorAxis, squeeze(RDiffAxes(i, :, 2))) / norm(curVigorAxis) / norm(squeeze(RDiffAxes(i, :, 2)));
            projVigorAxis2 = projVigorAxis * 10;
            scatter(curNeuralData_onPC(1), curNeuralData_onPC(2), 75, rewColors(k, :), 'filled', 'MarkerEdgeColor', diffColors(j, :), 'LineWidth', 1.5);
            plot([curNeuralData_onPC(1), curNeuralData_onPC(1) + projVigorAxis2(1)], ...
                            [curNeuralData_onPC(2), curNeuralData_onPC(2) + projVigorAxis2(2)], 'Color', rewColors(k, :));
            angle = acosd(norm(projVigorAxis));
            text(curNeuralData_onPC(1) + projVigorAxis2(1), curNeuralData_onPC(2) + projVigorAxis2(2), ...
                            num2str(angle), 'Color', rewColors(k, :));
        end
    end
%     curVigorAxis = vigorAxes_eachDir(i, :);
%     curVigorAxis = squeeze(curVigorAxis);
%     projVigorAxis(1) = dot(curVigorAxis, RDiffAxes(i, :, 1)) / norm(curVigorAxis) / norm(RDiffAxes(i, :, 1));
%     projVigorAxis(2) = dot(curVigorAxis, RDiffAxes(i, :, 2)) / norm(curVigorAxis) / norm(RDiffAxes(i, :, 2));
%     projVigorAxis2 = projVigorAxis * 10;
%     xlim = get(gca, 'XLim');
%     ylim = get(gca, 'YLim');
%     plot([xlim(1), xlim(1) + projVigorAxis2(1)], [ylim(1), ylim(1) + projVigorAxis2(2)], 'Color', 'k');
%     angle = acosd(norm(projVigorAxis));
%     text(xlim(1) + projVigorAxis2(1), ylim(1) + projVigorAxis2(2), num2str(angle), 'Color', 'k');
    hold off
    xlabel("Reward&Difficulty PC1 (spikes/s)");
    ylabel("Reward&Difficulty PC2 (spikes/s)");
    saveas(gcf, outputFolder + "vigorAxis_for_direction_" + directions(i) + ".png")
    close all;


    curVigorAxis = vigorAxes_eachDir(i, :);
    curVigorAxis = squeeze(curVigorAxis);
    RDiffAxis = squeeze(RDiffAxes(i, :, 2));
    angles(i, 1) = acosd(abs(dot(curVigorAxis, RDiffAxis) ./ norm(curVigorAxis) ./ norm(RDiffAxis)));
    angles(i, 2) = acosd(abs(dot(curVigorAxis, wDif(:, 1)) ./ norm(curVigorAxis) ./ norm(wDif(:, 1))));
end

% Allangles(d, :, :) = angles;
% end

% angles = Allangles(:, :, 1);
angles = angles(:, 1);

% plot the PDF of angles and PDF of null distribution
[nullAngles]  = randomAngle(19, 10000);

[counts, edges] = histcounts(nullAngles, 0:90, 'Normalization', 'pdf');
centers = (edges(1:end-1) + edges(2:end))/2;
figure; hold on;
bar(centers, counts, 'FaceColor', [0.7 0.7 0.7]);  % change color as needed
[counts, edges] = histcounts(angles, 0:90, 'Normalization', 'pdf');
centers = (edges(1:end-1) + edges(2:end))/2;
bar(centers, counts, 'FaceColor', [1 0 0]);  % change color as needed
xlabel('Angle');
ylabel('PDF');
title('Histogram with PDF on Y-axis');
saveas(gcf, outputFolderAngle + "RDiffvsNull.png");