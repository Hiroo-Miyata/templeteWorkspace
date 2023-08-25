clear all;
%% Input 
% dates = ["20220407"];
dates = ["20220216", "20220217", "20220218", "20220221", "20220222", "20220223", "20220225", "20220228", "20220301"];
for d = 1:length(dates)
date = dates(d);
rootDir = "../";
target = "HT";
beforeT = 200;
afterT = 375;
conditionState = [3 4 5 6];
minDelayTime = 100; 
outputFile = rootDir + "data/processed/non-stitched/" + date + "_"+target+"_"+num2str(beforeT)...
    +"_"+num2str(afterT)+"_"+"3456";

%% load
preprocessedFolder = rootDir + "data/preprocessed/" + date;
synchfile = dir(preprocessedFolder+'/*_TSK_*.mat');
neurofile = dir(preprocessedFolder+'/*_NER_*.mat');
% emgfile = dir(preprocessedFolder+'/*_EMG_*.mat');
% pupilfile = dir(preprocessedFolder+'/*_EYE_*.mat');
kinematicsfile = dir(preprocessedFolder+'/*_updatedKIN_*.mat');
load(synchfile.folder +"/"+ synchfile.name);
load(neurofile.folder +"/"+ neurofile.name);
% load(emgfile.folder +"/"+ emgfile.name);
% load(pupilfile.folder +"/"+ pupilfile.name);
load(kinematicsfile.folder +"/"+ kinematicsfile.name);

%% data around GoCue -200~+600 ms
if target == "TO"
    state = 3;
elseif target == "GC"
    state = 4;
elseif target == "HT"
    state = 6;
end

% array indicate extracted data
validTrials = false(length(trialData), 1);
%% extract data corresponding to the inputs
for i = 1:length(trialData)
    stateTransition = trialData(i).stateTable;
    trialNum = trialData(i).trial;
    %% if the trial has the state of interest: [3]
    if all(ismember(conditionState, stateTransition(1,:))) == 1 && trialNum > 24
        % if the trial is not bad trials
        if spikeInfo.badTrials(i) == 0
            TOidx = find(stateTransition(1, :)==3);
            % if the delay period is longer than 350 ms
            if minDelayTime < (stateTransition(2, TOidx+1) - stateTransition(2, TOidx))
                stateTime = stateTransition(2, find(stateTransition(1, :)==state));
                stateTime2 = find(trialData(i).time==stateTime);

                %gaussian kernel
                firingRate = zeros(size(neuralData(i).spikeMatrix)); % nneurons * times
                spikeTrain = neuralData(i).spikeMatrix;
%                 sigma = 25;
%                 kernel = normpdf(-3*sigma:3*sigma,0,sigma);
%                 for t=(1:size(firingRate, 1))
%                     firingRate(t, :) = 1000*conv(spikeTrain(t, :), kernel, "same");
%                 end
                firingRate = 1000 * spikeTrain;
                trialData(i).firingRates = firingRate(:, stateTime2-beforeT:stateTime2+afterT);
%                 trialData(i).emg = EMGData(i).signal(stateTime2-beforeT:stateTime2+afterT, :);
%                 trialData(i).goodEMGData = EMGData(i).goodEMGData;
                trialData(i).handKinematics = kinematicData(i);
                % GCTime = stateTransition(2, find(stateTransition(1, :)==4));
                % GCTime = find(trialData(i).time==GCTime);
                % TOTime = stateTransition(2, find(stateTransition(1, :)==3));
                % TOTime = find(trialData(i).time==TOTime);
                % trialData(i).pupilSize = pupilData(i).size(stateTime2-beforeT:stateTime2+afterT);
                % trialData(i).goodPupil = pupilData(i).goodPupil(stateTime2-beforeT:stateTime2+afterT);

                validTrials(i) = true;
            end
        end
    end
end

trialData = trialData(validTrials);
trialData = rmfield(trialData,["taskSynchTrialTime"]);
save(outputFile, "trialData")

end