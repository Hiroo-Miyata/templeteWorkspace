clear all;
%% Input 
dates = ["20220216", "20220217", "20220218","20220221","20220222","20220223","20220225","20220228","20220301"];

for d = 1:length(dates)
date = dates(d);
rootDir = "../";
target = "TO"; 
beforeT = 200;
afterT = 600;
conditionState = [3];
minDelayTime = 100; 
outputFile = rootDir + "data/processed/whole-stitched/" + date + "_"+target+"_"+num2str(beforeT)...
    +"_"+num2str(afterT)+"_"+"3";

%% load
preprocessedFolder = rootDir + "data/preprocessed/" + date;
synchfile = dir(preprocessedFolder+'/*_TSK_*.mat');
neurofile = dir(preprocessedFolder+'/*_StitchedWholeNER_*.mat');
spikeInfoFile = dir(preprocessedFolder+'/*_NER_*.mat');
emgfile = dir(preprocessedFolder+ '/*_EMG_*.mat'); 
% pupilfile = dir(preprocessedFolder+'/*_EYE_*.mat');
kinematicsfile = dir(preprocessedFolder+'/*_KIN_*.mat');
load(synchfile.folder +"/"+ synchfile.name);
load(neurofile.folder +"/"+ neurofile.name);
load(emgfile.folder +"/"+ emgfile.name);
load(spikeInfoFile.folder +"/"+ spikeInfoFile.name, "spikeInfo");
% load(pupilfile.folder +"/"+ pupilfile.name);
load(kinematicsfile.folder +"/"+ kinematicsfile.name);

%% data around GoCue -200~+600 ms
if target == "TO"
    state = 3;
elseif target == "GC"
    state = 4;
elseif target == "HT"
    state = 6;
elseif target == "SU"
    state = 7;
elseif target == "DF"
    state = 11;
elseif target == "RF"
    state = 12;
end

% array indicate extracted data
validTrials = false(length(trialData), 1);
%% extract data corresponding to the inputs
for i = 1:(length(trialData)-1)
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
                lastT = max(trialData(i).time);
                
                if (stateTime + afterT) < lastT
                    stateTime2 = find(trialData(i).time==stateTime);
                    %gaussian kernel
                    firingRate = zeros(size(neuralData(i).FScoreMatrix)); % nneurons * times
                    spikeTrain = neuralData(i).FScoreMatrix;
                    firingRate = spikeTrain;
                    trialData(i).firingRates = firingRate(:, stateTime2-beforeT:stateTime2+afterT);
                    
                    % GCTime = stateTransition(2, find(stateTransition(1, :)==4));
                    % GCTime = find(trialData(i).time==GCTime);
                    % TOTime = stateTransition(2, find(stateTransition(1, :)==3));
                    % TOTime = find(trialData(i).time==TOTime);
                    % trialData(i).pupilSize = pupilData(i).size(stateTime2-beforeT:stateTime2+afterT);
                    % trialData(i).goodPupil = pupilData(i).goodPupil(stateTime2-beforeT:stateTime2+afterT);
                    trialData(i).emg = EMGData(i).signal(stateTime2-beforeT:stateTime2+afterT, :);
                    trialData(i).goodEMGData = EMGData(i).goodEMGData;
                    
                    % there are 7 column. position, velocity, acceleration, speed, rotatedPosition, distanceFromCenter, distanceFromEndTarget
                    trialData(i).handKinematics.position = kinematicData(i).position(stateTime2-beforeT:stateTime2+afterT, 1:2);
                    trialData(i).handKinematics.velocity = kinematicData(i).velocity(stateTime2-beforeT:stateTime2+afterT, 1:2);
                    trialData(i).handKinematics.acceleration = kinematicData(i).acceleration(stateTime2-beforeT:stateTime2+afterT, 1:2);
                  

                    validTrials(i) = true;
                else
                    if trialData(i+1).trial == (trialNum + 1)
                        firingRate = zeros(size(neuralData(i).FScoreMatrix)); % nneurons * times
                        spikeTrain = neuralData(i).FScoreMatrix;
                        firingRate = spikeTrain;
                        [nneurons, ntimeBins] = size(firingRate);
                        combinedFR = nan(nneurons, afterT+beforeT+1);
                        dataLimit = lastT-stateTime;
                        combinedFR(:, 1:(dataLimit+beforeT+1)) = firingRate(:, end-dataLimit-beforeT:end);
                        firingRate = zeros(size(neuralData(i+1).FScoreMatrix)); % nneurons * times
                        spikeTrain = neuralData(i).FScoreMatrix; firingRate = spikeTrain;
                        combinedFR(:, (dataLimit+beforeT+2):end) = firingRate(:, 1:afterT-dataLimit);
                        trialData(i).firingRates = combinedFR;

                        combineEMG = nan(afterT+beforeT+1, 5);
                        combineEMG(1:(dataLimit+beforeT+1), :) = EMGData(i).signal(end-dataLimit-beforeT:end, :);
                        combineEMG((dataLimit+beforeT+2):end, :) = EMGData(i+1).signal(1:afterT-dataLimit, :);
                        trialData(i).emg = combineEMG;
                        trialData(i).goodEMGData = EMGData(i).goodEMGData & EMGData(i+1).goodEMGData;

                        trialData(i).handKinematics.position = nan(afterT+beforeT+1, 2);
                        trialData(i).handKinematics.position(1:(dataLimit+beforeT+1), :) = kinematicData(i).position(end-dataLimit-beforeT:end, 1:2);
                        trialData(i).handKinematics.position((dataLimit+beforeT+2):end, :) = kinematicData(i+1).position(1:afterT-dataLimit, 1:2);
                        trialData(i).handKinematics.velocity = nan(afterT+beforeT+1, 2);
                        trialData(i).handKinematics.velocity(1:(dataLimit+beforeT+1), :) = kinematicData(i).velocity(end-dataLimit-beforeT:end, 1:2);
                        trialData(i).handKinematics.velocity((dataLimit+beforeT+2):end, :) = kinematicData(i+1).velocity(1:afterT-dataLimit, 1:2);
                        trialData(i).handKinematics.acceleration = nan(afterT+beforeT+1, 2);
                        trialData(i).handKinematics.acceleration(1:(dataLimit+beforeT+1), :) = kinematicData(i).acceleration(end-dataLimit-beforeT:end, 1:2);
                        trialData(i).handKinematics.acceleration((dataLimit+beforeT+2):end, :) = kinematicData(i+1).acceleration(1:afterT-dataLimit, 1:2);
                        


                        validTrials(i) = true;
                    end
                end
%                 trialData(i).emg = EMGData(i).signal(stateTime2-beforeT:stateTime2+afterT, :);
%                 trialData(i).goodEMGData = EMGData(i).goodEMGData;
%                 trialData(i).handKinematics = kinematicData(i);
            end
        end
    end
end

trialData = trialData(validTrials);
trialData = rmfield(trialData,["taskSynchTrialTime"]);
save(outputFile, "trialData");

end