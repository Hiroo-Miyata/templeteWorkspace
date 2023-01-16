function [neuralData, spikeInfo] = p1_spikePreprocessing(trialData, neuralData,minSpikeRate,corrThresh, plotoutputFolder)
% This function is the 1st processing step after spike-sorted data from the
% Batista lab has been run through initial preprocessing. Here, we do three
% things:
%
% (1) Get full trial spike rates (delay through end of reach)
% (2) Remove outlier units/trials; like HARD outlier, like "wtf something
% is wrong with the data here" outlier. Also if it's the wrong task
% (3) Remove crosstalk units (i.e. appeared duplicate unit)
%
% Inputs:
% - trialData: struct.array task information of each trial 
% - neuralData: cell.array neural spike data of each trial 
% - minSpikeRate: scalar in Hz that represents the minimum average firing
% rate of a neuron across all trials to be worth keeping. We use this for
% minimum averaging firing rate across neurons for individual trials as
% well.
% - corrThresh: scalar that is the maximum correlation allowed between two
% units' spike trains before they're considered one in the same (and hence,
% one must be removed).
%
%
% Outputs:
% - spikeInfo struct include badTrials and badUnits information
% Hiroo Miyata, 11/21/2022

%% info of Rocky dataset
targetOnsetState = 3;
endStates = [7 10 11 12 13];
nunits = size(neuralData(1).spikeMatrix,1);
ntrials = length(neuralData);


%% (1) Get full trial spike rates
badTrials = false(ntrials,1);
badUnits = false(nunits,1);
fullTrialSpikeRateMat = nan(nunits,ntrials);
for j = 1:ntrials
    trans = trialData(j).stateTable;
    startTime = trans(2,find(trans(1,:)==targetOnsetState,1));
    endTime = min(trans(2,find(ismember(trans(1,:),endStates),1)), ...
        max(trialData(j).time)); % For S trials, sometimes a few ms got chopped off, so "reward" time is beyond neural data
    if isempty(startTime) || isempty(endTime)
        badTrials(j) = true;
    else % We assume neural data is on 1 ms basis
        fullTrialSpikeRateMat(:,j) = 1000*mean(neuralData(j).spikeMatrix(:, startTime:endTime),2);
    end
end; clear j

firingRateAll = fullTrialSpikeRateMat; % tmp file for plot
    
%% First, if they don't spike enough, they're out! This probs means the
% trial is an artifact of some sort
noSpikeTrials = mean(fullTrialSpikeRateMat) < minSpikeRate;
badTrials(noSpikeTrials) = true;

% One artifact that also seems to happen is like 5 neurons have
% ridiculously high FR, the rest are 0 - so the trial avg over neurons
% looks normal. We'll just see if the median firing rate for a trial is
% less than the minimum spiking rate, indicating silence on > half of
% the units.
sparseNeuronTrials = median(fullTrialSpikeRateMat) < minSpikeRate;
badTrials(sparseNeuronTrials) = true;

% Next, we'll identify and remove "outlier trials". To do this, we 
% go through each trial individually, calculate the median absolute 
% deviation (MAD) of the other trials for each neuron, and see if the 
% current trial is outside of +/- 3 MADs for >= 75% of the neurons. 
% If so, we call it an outlier. This may re-pick up the no-spike trials
% above - that's fine.
outlierTrials = false(ntrials,1);
noutlierUnits = nan(ntrials,1);
for j = 1:ntrials
    if ~badTrials(j)
        frOfTrial = fullTrialSpikeRateMat(:,j);
        indsToNotUse = badTrials;       % which inds won't be in the MAD calc?
        indsToNotUse(j) = true;         % we also shouldn't use the current ind
        otherData = fullTrialSpikeRateMat(:,~indsToNotUse);
        otherDataMedian = median(otherData,2);
        curMADs = mad(otherData,1,2);   % 1 = median, 2 = dimension
        skipUnits = curMADs==0;         % skip these units; they'll probably be removed later
        maxDeviation = curMADs(~skipUnits)*3;
        outlierUnits = abs(otherDataMedian(~skipUnits)-frOfTrial(~skipUnits)) >= maxDeviation;
        noutlierUnits(j) = sum(outlierUnits);
        if sum(outlierUnits)>= 0.75*sum(~skipUnits)
            outlierTrials(j) = true;
        end
    end
end; clear j

% Another way - use the euclidean distance from the median; only used
% for Rocky so far
trialDistFromMedian = sqrt(sum((fullTrialSpikeRateMat-median(fullTrialSpikeRateMat,2)).^2));
distMedian = median(trialDistFromMedian);
distMAD = mad(trialDistFromMedian,1);
outlierTrials = outlierTrials | (abs(trialDistFromMedian-distMedian) > 4*1.48*distMAD)'; % this is roughly equivalent to 4 stds away if it was actually normal

% Just for counting, get reward labels and # of jackpots removed
rewardLabels = [trialData.rewardLabel];


% Remove all bad trials
badTrials(outlierTrials) = true;
disp(['# trials removed: ' num2str(sum(badTrials))])
disp(['# of Jackpots removed: ' num2str(sum(rewardLabels(badTrials)==4))])

fullTrialSpikeRateMat(:,badTrials) = [];    % remove it from these too
spikeInfo.badTrials = badTrials;            % removes bad trials...nice and easy.
spikeInfo.metadata = zeros(length(neuralData(1).channel),8);
spikeInfo.metadata(:,1) = neuralData(1).channel;
spikeInfo.metadata(:,2) = neuralData(1).sort;


%% (2b) Remove outlier and unsorted units
disp('Removing bad units')
% First, if they don't spike enough, they're out! They'll probs break
% the FA model...
noSpikeUnits = mean(fullTrialSpikeRateMat,2) < minSpikeRate;
badUnits(noSpikeUnits) = true;
spikeInfo.metadata(:,3) = badUnits;

% This ^ captures units who fire too slowly for the FA model; however,
% the other prime concern is "units" that are just artifacts -
% typically, some ridiculously high FR with virtually 0 variance. We
% don't want these units in our data! To avoid them, we'll take the
% stdev of the data; if a unit is < 1 Hz across all trials, remove it.
lowVarUnits = std(fullTrialSpikeRateMat,[],2) < minSpikeRate;
badUnits(lowVarUnits) = true;
spikeInfo.metadata(:,4) = badUnits;

% Occasionally a unit will be present / active for like half
% of the session, then disappear. This type of artifact can really
% throw off dimensionality reduction, so let's clean it up. To check
% this, we'll simply look at the mean of the first vs. last 1/4 of the
% session's trials.
quarterOfSession = round(size(fullTrialSpikeRateMat, 2)/4);
earlyFR = mean(fullTrialSpikeRateMat(:,1:quarterOfSession),2);
lateFR = mean(fullTrialSpikeRateMat(:,end-quarterOfSession+1:end),2);
droppedUnits = (earlyFR-lateFR)./earlyFR > 0.75; % >75% drop in FR
badUnits(droppedUnits) = true;
spikeInfo.metadata(:,5) = badUnits;

% We can also check for the opposite
lateUnits = (lateFR-earlyFR)./lateFR > 0.75; % >75% increase in FR (but not for super low FR units)
badUnits(lateUnits) = true;
spikeInfo.metadata(:,6) = badUnits;

% Another artifact is that randomly the channel just goes
% crazy and gets ridiculously high values. If this happens, it seems to
% hold for a while, so we'll look for some sustained ridiculously high
% activity.
for i=(1:nunits)
    smoothedTrialFR(i, :) = conv(fullTrialSpikeRateMat(i, :),1/50*ones(50,1),'valid'); % running 50 trial window
end
artifactHighUnits = any(smoothedTrialFR > 200,2); % if the avg of 50 trials is ever > 200Hz, there's a problem!
badUnits(artifactHighUnits) = true;
spikeInfo.metadata(:,7) = badUnits;

%% I couldn't understand. Finally, identify sort ID, and if it's not between 1-6, throw it out.
% unsortedUnits = ~ismember(data(1).TrialData.neuralData.sort,1:6);
% badUnits(unsortedUnits) = true;

%% (3) Remove duplicated "cross-talk" units
binSize = 10;
dataForCrosstalk = [];
badTrialRemovedNeuralData = neuralData;
badTrialRemovedNeuralData(badTrials) = [];
for j = 1:min(20,length(badTrialRemovedNeuralData)) % 10 trials should be enough
    curTrial = badTrialRemovedNeuralData(j).spikeMatrix;
    curTrialBinned = [];
    curInd = 1;
    while curInd+binSize <= size(curTrial,1)
        curTrialBinned(:,end+1) = sum(curTrial(:, curInd:curInd+binSize),2);
        curInd = curInd+binSize;
    end
    dataForCrosstalk = [dataForCrosstalk curTrialBinned];
end

corrMat = corr(dataForCrosstalk');
corrMatNoDiag = corrMat-diag(diag(corrMat));
corrMatNoDiag(badUnits,:) = 0; % used to exclude already-labeled bad units from this analysis
corrMatNoDiag(:,badUnits) = 0;
neuronCrosstalkSums = sum(corrMatNoDiag > corrThresh);
crosstalkUnits = [];
while sum(neuronCrosstalkSums)~=0
    [~,worstOffender] = max(neuronCrosstalkSums); % which neuron has the most crosstalks? Kill it.
    crosstalkUnits(end+1) = worstOffender;
    corrMatNoDiag(worstOffender,:) = 0;
    corrMatNoDiag(:,worstOffender) = 0;
    neuronCrosstalkSums = sum(corrMatNoDiag > corrThresh);
end
badUnits(crosstalkUnits) = true;
spikeInfo.metadata(:,8) = badUnits;

% Remove all bad units
disp(['# units removed: ' num2str(sum(badUnits))])
spikeInfo.badUnits = badUnits;
spikeInfo.metadata = array2table(spikeInfo.metadata, ...
    'VariableNames', {'Ch', 'Unit', 'minSpike', 'minSpikeVar', ...
    '75%drop', '75%increase', 'crazyChannel', 'crossTalk'});
for j = 1:ntrials
    neuralData(j).channel(badUnits) = [];
    neuralData(j).sort(badUnits) = [];
    neuralData(j).spikeMatrix(badUnits, :) = [];
    neuralData(j).meanWaveforms(badUnits, :) = [];
end; clear j;

plotHeatMap(firingRateAll, badTrials, badUnits, plotoutputFolder)
disp('p1_spikePreprocessing complete!')
end

