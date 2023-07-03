close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];
trialDataAll = struct.empty;
for date = dates
    load("../data/processed/stitched/2022"+date+"_GC_200_600_3456.mat");
    taskfile = dir("../data/preprocessed/2022"+date+ "/*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    trialNum = [trialData.trial];

    behavData = behavData(ismember([behavData.trial], trialNum));
    for i = 1:length(trialData)
        trialData(i).reactionTime = behavData(i).t4_reactionTime - behavData(i).t3_goCueTime;
        % add delay time
        trialData(i).delayTime = behavData(i).t3_goCueTime - behavData(i).t2_targetOnsetTime;
    end
    trialDataAll = cat(1, trialDataAll, trialData);
end
trialData = trialDataAll;
outputFolderBehav = "../results/202306w3-results-summary-all/neuralData_GoCue/motorPlanning_targetPlane/";

% ndates = length(dates);
% for d = 1:ndates
% date = dates(d);
% load("../data/processed/not-stitched/2022"+date+"_GC_200_600_3456.mat");
% outputFolderBehav = "../results/202306w3-results-summary-all/neuralData_GoCue/motorPlanning_targetPlane/"+date;

%% get only the trials with shorter delay time < 650ms
% delayTimes = [trialData.delayTime];
% trialData = trialData(delayTimes < 650);


% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};

% Get the average firing rates of each neuron in each reward and difficulty condition
analysisBin = (50:250); % Pay Attention!! GC=(50:250), TO=(100:300)
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


meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[2 3]));
[wTP,zTP,zTPeigVls_TP,~,~] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndirections-1);
save("../interim/GC_targetPlaneAxis_all.mat", "wTP", "zTP", "zTPeigVls_TP");

neuralData_onDPCs = (neuralData - mean(neuralData, 1)) * wTP;

%% signal variance and noise variance
tpSize = zeros(nrewards, ndifficulties);
covSize = zeros(nrewards, ndifficulties);
noiseVariances = zeros(nrewards, ndifficulties);
for rew=1:nrewards
    for diff=1:ndifficulties
        x=zeros(length(directions), 1);
        y=zeros(length(directions), 1);
        nV = zeros(ndirections, 1);
        for dire = 1:ndirections
            curInds = rewardLabels==rewards(rew) & difficultyLabels==difficulties(diff) & directionLabels==directions(dire);
            tmp = neuralData_onDPCs(curInds, 1:2);
            x(dire) = mean(tmp(:, 1));
            y(dire) = mean(tmp(:, 2));
            nV(dire) = trace(cov(tmp));
        end
        x = x - mean(x);
        y = y - mean(y);
        tpSize(rew, diff) = mean(rssq([x y], 2));
        covSize(rew, diff) = det(cov(x, y));
        noiseVariances(rew, diff) = mean(nV);
    end
end

%% plotting
VS_reward_size_calculated(covSize, "Label", "Signal Variance in Target Plane (a.u.)", "OutputFolder", outputFolderBehav+"SignalVariance")
VS_reward_size_calculated(noiseVariances, "Label", "Noise Variance in Target Plane (a.u.)", "OutputFolder", outputFolderBehav+"NoiseVariance")
