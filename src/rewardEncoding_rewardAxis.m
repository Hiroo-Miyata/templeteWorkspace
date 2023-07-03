close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
ndates = length(dates);
Yall = cell(ndates, 3, 2, 8);
for d = 1:ndates
date = dates(d);
rootDir = "../";
load(rootDir+"data/processed/non-stitched/2022"+date+"_GC_200_600_3456.mat");
% load(rootDir+"data/processed/stitched-whole/2022"+date+"_TO_200_700_3.mat");
% load(rootDir+"data/processed/non-stitched/2022"+date+"_HT_200_375_3456.mat");
outputFolder = "../results/202306w3-results-summary-all/neuralData_GoCue/projection-2D/"; makeDir(outputFolder); 
% outputFolderProjection = "../results/202306w3-results-summary-all/neuralData_GoCue/rewardAxis-projection/";
% outputFolderStretching = "../results/202306w3-results-summary-all/neuralData_GoCue/stretching/";
analysisBin = (50:250); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(400:600)
axisOutputFileName = rootDir+"/interim/nonstitched_GC_rewardAxis_"+date+".mat";

% remove calibration part 
trialNums = [trialData.trial];
trialData = trialData(trialNums > 24);
% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trial]; 
trialStatusLabels = [trialData.trialStatus];

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

rewardAxisName = "rewardAxisBoth";
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
for i=1:nrewards
    for j=1:ndifficulties
        for k=1:ndirections
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
            Yall{d, i, j, k} = neuralData_onDPC(curInds, 1);
        end
    end
end

%% angle of wDif and wR
dot_product = dot(wRTiny(:,1), wRHuge(:,1));
angle = acos(dot_product) * 180 / pi;
if angle > 90
    angle = 180 - angle;
end
disp(angle);

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
% maxlim = ceil(max(X(:)) / 10)*10;
% minlim = floor(min(X(:)) / 10)*10;
% ylim([minlim maxlim]); xlim([minlim maxlim]);
ylim([-30 20]); xlim([-30 20]);
xlabel("RDiff Axis: PC1 (Hz): "+round(eigVls_RDiff(1)/sum(eigVls_RDiff)*100, 2)+"%");
ylabel("RDiff Axis: PC2 (Hz): "+round(eigVls_RDiff(2)/sum(eigVls_RDiff)*100, 2)+"%");
legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
saveas(gcf, outputFolder+"2Dprojection-"+date+".jpg");


%% plot reward PC1 and difficulty PC1
% rAxes = ["Tiny", "Huge", "All"];
% for k=1:3 % 1:2
%     Y = zeros(3,2); % reward x difficulty
%     Yerr = zeros(3,2);
%     if k==1
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRTiny(:,1);
%     elseif k==2
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRHuge(:,1);
%     else
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
%     end
% 
%     for i=1:nrewards
%         for j=1:ndifficulties
%             curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j);
%             Y(i,j) = mean(neuralData_onDPC(curInds));
%             Yerr(i,j) = std(neuralData_onDPC(curInds))/sqrt(sum(curInds));
%         end
%     end
%     figure;
%     for j=1:ndifficulties
%         errorbar(1:nrewards, Y(:,j), Yerr(:,j), "Color", diffColors(j,:), 'LineWidth', 2); hold on;
%     end
%     set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
%     xticks([1 2 3]); xticklabels(["S", "M", "L"]); xlim([0.7 3.3]); legend(["Tiny", "Huge"], Location="best")
%     set(gcf,'position',[100,100,300,650]);
%     ylabel("Reward Axis (Hz)");
%     saveas(gcf, outputFolderProjection+date+"-Axis-in-"+rAxes(k)+".jpg"); close all;
% end


%% get the gap between the two difficulty conditions in each reward condition
% Yall(d, r) = projection of the mean firing rate in the Tiny difficulty condition on the reward axis - that in the Huge difficulty condition

% Y = zeros(3,2); % reward x difficulty
% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
% for i=1:nrewards
%     for j=1:ndifficulties
%         curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & trialStatusLabels == 1;
%         Y(i,j) = mean(neuralData_onDPC(curInds));
%     end
% end
% Yall(d, :) = Y(:,1) - Y(:,2);

%% the reward axis value in delay failed trials

% prepare matrix for reward axis value in delay failed trials and non-delay failed trials
% rewardAxisValue = zeros(3, 2 ,2); % reward x difficulty x delay failed/non-delay failed
% rewardAxisValueErr = zeros(3, 2 ,2); % reward x difficulty x delay failed/non-delay failed
% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
% 
% delayFailedLabels = zeros(size(rewardLabels));
% for i=1:length(trialData)
%     stateTable = trialData(i).stateTable;
%     if all(ismember([3 11], stateTable(1,:))) == 1
%         delayFailedLabels(i) = 1;
%     end
% end
% disp(sum(delayFailedLabels)/length(delayFailedLabels));
% 
% for i=1:nrewards
%     for j=1:ndifficulties
%         for k=1:2 % non-delay / delay
%             curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & delayFailedLabels == k-1;
%             rewardAxisValue(i,j,k) = mean(neuralData_onDPC(curInds));
%             rewardAxisValueErr(i,j,k) = std(neuralData_onDPC(curInds))/sqrt(sum(curInds));
%         end
%     end
% end
% figure;
% for j=1:ndifficulties
%     h(1) = plot(1:nrewards, rewardAxisValue(:,j,1), "Color", diffColors(j,:), 'LineWidth', 2); hold on;
%     h(2) = plot(1:nrewards, rewardAxisValue(:,j,2), "Color", diffColors(j,:), 'LineWidth', 2, 'LineStyle', '--'); hold on;
% end
% 
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% xticks([1 2 3])
% xticklabels(["S", "M", "L"])
% xlim([0.7 3.3])
% legend(h, ["non-delay failed", "delay failed"], Location="best")
% set(gcf,'position',[100,100,300,650]);
% ylabel("Reward Axis (Hz)");
% saveas(gcf, outputFolder+"Effect-DelayFailure-"+date+".jpg");

%% time trajectory 
% figure;
% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
% for r=1:nrewards
%     for diff=1:ndifficulties
%         curInds = rewardLabels == rewards(r) & difficultyLabels == difficulties(diff);
%         X = trialNums(curInds);
%         Y = movmean(neuralData_onDPC(curInds), 50);
%         p(diff) = plot(X, Y , Color=diffColors(diff,:), LineWidth=2);hold on;
%         scatter(X([1 end]), Y([1 end]), 100, rewColors(r, :), "filled"); hold on;
%     end
% end
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% legend(p, ["Tiny", "Huge"], Location="best")
% xlabel("Trials")
% ylabel("Reward Axis: PC1 (Hz)");
% saveas(gcf, outputFolder+"trajectory-"+date+".jpg");

%% plot PSTH projected on reward axis
% ntimeBin = size(neuralActivity, 2);
% PSTHs_onPCs = zeros(nrewards, ndifficulties, ntimeBin);

% for i = 1:nrewards
%     for j = 1:ndifficulties
%         curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
%         PSTH = squeeze(mean(neuralActivity(:, :, curInds), 3));
%         sigma = 25;
%         kernel = normpdf(-3*sigma:3*sigma,0,sigma);
%         for t=(1:size(PSTH, 1))
%             PSTH(t, :) = conv(PSTH(t, :), kernel, "same");
%         end
%         PSTH = PSTH';
%         PSTHs_onPCs(i, j, :) = PSTH * wR(:,1);
%     end
% end

% figure;
% hold on;
% for i = 1:nrewards
%     for j = 1:ndifficulties
%         PSTH = squeeze(PSTHs_onPCs(i, j, :));
%         plot(1:ntimeBin, PSTH, Color=diffColors(j,:), LineWidth=2); hold on;
%         scatter([1 ntimeBin], [PSTH(1) PSTH(ntimeBin)], 100, rewColors(i, :), 'filled');
%     end
% end
% xlabel("Time (ms)");
% ylabel("Projection on reward axis");
% xticks([1 201 ntimeBin]);
% xticklabels(["-200", "GC", num2str(ntimeBin - 201)]);
% set(gcf,'position',[100,100,450,650]);
% saveas(gcf, outputFolder+"Trajectory-rewardaxis-"+date+".jpg");

%% angle of wR at TO and wR at GC (same script is applied to different file)
% load(rootDir+"data/processed/"+date+"_GC_200_600_3456.mat");
% outputFolder = "../results/202301w1/rewardAxisTrajectoryAtTO/"+date; %rewardAxisTrajectoryAtTO
% directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
% rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
% difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
% rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
% diffColors = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
% analysisBin = (50:250);
% neuralActivity = cat(3, trialData(:).firingRates);
% neuralData = mean(neuralActivity(:, analysisBin, :), 2);
% neuralData = squeeze(neuralData)';
% [ntrials,nneurons] = size(neuralData);
% neuralData_byParameters_mean = nan(ndirections,nrewards,nneurons);
% for d = 1:ndirections
%     for r = 1:nrewards
%         curInds = directionLabels==directions(d) & rewardLabels==rewards(r);
%         neuralData_byParameters_mean(d,r,:) = mean(neuralData(curInds,:));
%     end; clear r
% end; clear d
% rewMeans = squeeze(mean(neuralData_byParameters_mean,1)); % nrewards x nneurons
% muR = mean(rewMeans);
% [wR2,zR2,eigVls_R2,~,pve_R2] = pca(rewMeans-muR,'numcomponents',nrewards-1);
% dot_product = dot(wR(:,1), wR2(:,1));
% angle = acos(dot_product) * 180 / pi;
% disp(angle);


% FR_mean = zeros(nneurons, nrewards);
% FR_std = zeros(nneurons, nrewards);
% FRs = squeeze(mean(neuralActivity(:, analysisBin, :), 2)); % nneuron x ntrial
% for r = 1:nrewards
%     FR_mean(:,r) = mean(FRs(:,rewardLabels==rewards(r)),2);
%     FR_std(:,r) = std(FRs(:,rewardLabels==rewards(r)),0,2);
% end; clear r
% 
% rewardNeuronsCount = zeros(4,1);
% rewardNeurons = cell(4,1);
% 
% for n = 1:nneurons
%     dif_SM = FR_mean(n,1) - FR_mean(n,2);
%     dif_ML = FR_mean(n,2) - FR_mean(n,3);
%     if (dif_SM > 0 && dif_ML > 0) || (dif_SM < 0 && dif_ML < 0)
%         % get p-value the data size is different in each reward condition
%         [~,p_SM] = ttest2(FRs(n,rewardLabels==rewards(1)), FRs(n,rewardLabels==rewards(2)));
%         [~,p_ML] = ttest2(FRs(n,rewardLabels==rewards(2)), FRs(n,rewardLabels==rewards(3)));
%         if p_SM < 0.05 && p_ML < 0.05
%             rewardNeuronsCount(1) = rewardNeuronsCount(1) + 1;
%             rewardNeurons{1} = [rewardNeurons{1} n];
%         elseif p_SM < 0.05 && p_ML >= 0.05
%             rewardNeuronsCount(2) = rewardNeuronsCount(2) + 1;
%             rewardNeurons{2} = [rewardNeurons{2} n];
%         elseif p_SM >= 0.05 && p_ML < 0.05
%             rewardNeuronsCount(3) = rewardNeuronsCount(3) + 1;
%             rewardNeurons{3} = [rewardNeurons{3} n];
%         else
%             rewardNeuronsCount(4) = rewardNeuronsCount(4) + 1;
%             rewardNeurons{4} = [rewardNeurons{4} n];
%         end
%     end
% end; clear n
% 
% disp(rewardNeuronsCount)
% disp(nneurons)


close all;
end
rewardAxisRaw = Yall;
rewardAxisMean = zeros(8,3,2,8);
rewardAxisStd  = zeros(8,3,2,8);
for d=1:8
    for j=1:ndifficulties
        for i=1:nrewards
            for k=1:ndirections
                Ys = Yall{d,i,j,k};
                rewardAxisMean(d,i,j,k) = mean(Ys);
                rewardAxisStd(d,i,j,k)  = std(Ys) / sqrt(length(Ys));
            end
        end
    end
end


comprisons = [3 1; 3 2; 2 1];
comprisonsName = ["L_vs_S", "L_vs_M", "M_vs_S"];
for k = 1:3

    X = zeros(ndates,ndirections);
    Y = zeros(ndates,ndirections);
    for i = 1:ndates
        for j = 1:ndirections
            X(i,j) = rewardAxisMean(i,comprisons(k, 1),1,j) - rewardAxisMean(i,comprisons(k, 2),1,j);
            Y(i,j) = rewardAxisMean(i,comprisons(k, 1),2,j) - rewardAxisMean(i,comprisons(k, 2),2,j);
        end
    end
%     rankTest2d(X(:), Y(:), "XLabel", "Reward axis range, Tiny targets (Spikes/s)", ...
%                 "YLabel", "Reward axis range, Huge targets (Spikes/s)", ...
%                 "OutputFolder", outputFolderStretching+rewardAxisName+"-stretch-"+comprisonsName(k))
end