close all;clear all;
addpath(genpath("./utils/plotting/"));
addpath(genpath("./utils/processing/"));
% date = "0413";
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
Yall = zeros(length(dates), 2, 8, 2);
for d = 1:length(dates)
date = dates(d);
outputFolder = "../results/202306w3-results-summary-all/neuralData_GoCue/motorPlanning_afsharAxis/"; makeDir(outputFolder);

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
reactionTimes = [trialData.reactionTime];
curInds = ~isnan(reactionTimes) & reactionTimes > 200 & reactionTimes < 450;
trialData = trialData(curInds);

% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNum = [trialData.trial];
reactionTimes = [trialData.reactionTime];

% Get mean neural trajectory(p_t) of each direction (average across rewards and target sizes)
neuralActivity = cat(3, trialData(:).firingRates);
[neuralActivity, nneurons, ntimebins, ntrials] = kernelSmooth(neuralActivity);

% do the analysis by k-fold cross validation
kfold = 2;

Yalpha = zeros(kfold, ndirections, 2);
Ybeta  = zeros(kfold, ndirections, 2);
rTcorr = zeros(kfold, ndirections, 2);
for f = 1:kfold
    % get the training and testing data
    curInds = mod(1:ntrials, kfold) == f-1;
    testingData = neuralActivity(:, :, curInds);
    newRewL = rewardLabels(curInds);
    newDifL = difficultyLabels(curInds);
    newDirL = directionLabels(curInds);
    newRT = reactionTimes(curInds);
    % get the mean neural trajectory of each direction in training data
    meanNeuralTrajectory = zeros(ndirections, nneurons, ntimebins);
    for i = 1:ndirections
        curInds2 = directionLabels==directions(i) & ~curInds;
        meanNeuralTrajectory(i, :, :) = mean(neuralActivity(:, :, curInds2), 3);
    end
    % Get the vector p = p_GC+100ms – p_GC [nneurons * 1] in each direction
    % GC = 201; and 1 data point means 1 milisecond
    GCvector = zeros(ndirections, nneurons);
    for i = 1:ndirections
        GCvector(i, :) = squeeze(meanNeuralTrajectory(i, :, 201+100) - meanNeuralTrajectory(i, :, 201));
        GCvector(i, :) = GCvector(i, :) / norm(GCvector(i, :));
    end

    % for i = 1:ndirections-1
    %     for j = i+1:ndirections
    %         dot_product = dot(GCvector(i, :), GCvector(j, :));
    %         angle = acos(dot_product) * 180 / pi;
    %     end
    % end

    % Project single trial data (v) on p: α = dot((v_GC – p_GC),p)
    alpha = zeros(sum(curInds), 1);
    for i = 1:sum(curInds)
        alpha(i) = dot(squeeze(testingData(:, 201, i))' - squeeze(meanNeuralTrajectory(newDirL(i), :, 201)), GCvector(newDirL(i), :));
    end

    % project the vector of neural speed from Gc -10ms to Gc+10ms on the vector p
    beta = zeros(sum(curInds), 1);
    for i = 1:sum(curInds)
        beta(i) = dot(squeeze(testingData(:, 201+10, i))' - squeeze(testingData(:, 201-10, i))', GCvector(newDirL(i), :));
    end 

    %% plot 
    % alpha is z-indexed in each direction: Y
    % Plot the histogram of the alpha in each difficulty condition with different colors
    % color is diffColors

    for i = 1:ndirections
        % if d == 1 & f == 1
        %     % plot the distribution in each reward / each difficulty
        %     figure;
        %     for j = 1:ndifficulties
        %         curInds2 = newDifL==difficulties(j) & newDirL==directions(i);
        %         h(j) = histogram(alpha(curInds2), -6:0.2:6, 'FaceColor', diffColors(j, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        %         hold on;
        %         % add median line
        %         line([median(alpha(curInds2)) median(alpha(curInds2))], ylim, 'Color', diffColors(j, :), 'LineWidth', 2);
        %     end
        %     saveas(gcf, outputFolder+"zzz-alpha-distribution-TH.png");

        %     figure;
        %     for r = 1:nrewards
        %         curInds2 = newRewL==rewards(r) & newDirL==directions(i);
        %         h(r) = histogram(alpha(curInds2), -6:0.2:6, 'FaceColor', rewColors(r, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        %         hold on;
        %         % add median line
        %         line([median(alpha(curInds2)) median(alpha(curInds2))], ylim, 'Color', rewColors(r, :), 'LineWidth', 2);
        %     end
        %     saveas(gcf, outputFolder+"zzz-alpha-distribution-SML.png");
        %     close all;
        % end 
        curInds2 = newDirL==directions(i) & ~isnan(alpha') & ~isnan(beta');
        zalpha = zscore(alpha(curInds2));
        zbeta = zscore(beta(curInds2));
        
        Yalpha(f, i, 1) = mean(zalpha(newRewL(curInds2)==3)) - mean(zalpha(newRewL(curInds2)==1));
        Yalpha(f, i, 2) = mean(zalpha(newDifL(curInds2)==difficulties(2))) - mean(zalpha(newDifL(curInds2)==difficulties(1)));
        Ybeta(f, i, 1) = mean(zbeta(newRewL(curInds2)==3)) - mean(zbeta(newRewL(curInds2)==1));
        Ybeta(f, i, 2) = mean(zbeta(newDifL(curInds2)==difficulties(2))) - mean(zbeta(newDifL(curInds2)==difficulties(1)));

        rTcorr(f, i, 1) = corr(alpha(curInds2), newRT(curInds2).');
        rTcorr(f, i, 2) = corr(beta(curInds2), newRT(curInds2).');
    end
end
GCvector = GCvector.';
save("../interim/GC_AfsharAxis_"+date+".mat", "GCvector")
Yall(d, :, :, :) = rTcorr;

% for i = 1:2
%     if i == 1
%         Ys = Yalpha;
%         metricsStr = "position";
%     else
%         Ys = Ybeta;
%         metricsStr = "speed";
%     end
%     for s = 1:2
%         if s == 1
%             Y = Ys(:, :, 1);
%             titleStr = "Reward";
%         else
%             Y = Ys(:, :, 2);
%             titleStr = "Target Size";
%         end
% 
%         % corrHistogram(Y(:), "Label", titleStr+" Effect: " + metricsStr, "Title", titleStr+" vs AfsharAxis-at-GC", "OutputFolder", outputFolder+metricsStr+"-"+titleStr)
%     end
% end
close all;

end
Y = Yall(:, :, :, 1);
corrHistogram(Y(:), "Label", "Correlation Coefficient", "Title", "RT vs AfsharAxis-at-GC", "OutputFolder", outputFolder+"all");