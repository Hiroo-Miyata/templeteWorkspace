close all;clear all;
% date = "0413";
dates = ["0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"];

for d = 1:length(dates)
date = dates(d);
% outputFolder = "../results/202302w2/ReactionTime-DelayTime/";

%% get reaction time
dataFolder = "../data/preprocessed/2022"+date+"/";
taskfile = dir(dataFolder+ "*" +date+ "*_TSK_*.mat");
load(taskfile.folder +"/"+ taskfile.name);
kinfile = dir(dataFolder+ "*" +date+ "*_KIN_*.mat");
load(kinfile.folder +"/"+ kinfile.name);
trialData = arrayfun(@(x) setfield(trialData(x), "handKinematics", kinematicData(x)),1:length(trialData));
validTrialData = struct.empty(0);

for i=(1:length(trialData))
    stateTransition = trialData(i).stateTable;
    if all(ismember([3 4], stateTransition(1,:))) == 1
        validTrialData = cat(1, validTrialData, trialData(i));
    elseif all(ismember([3 11], stateTransition(1,:))) == 1
        validTrialData = cat(1, validTrialData, trialData(i));
    else
        
    end
end
tasktmp = taskInfo.taskParams;
tasktmp.subjectName = taskInfo.subjectName;
updatedTrialData = p5_updateBehavior(validTrialData, tasktmp);


behavData = struct.empty(0);
kinematicData = struct.empty(0);
trialNum = [updatedTrialData.trial];
% the length is the same as the number of trialData
% store updateTrialData if the trial is valid

for i = 1:length(trialData)
    kinematicData(i).trial = trialData(i).trial;
    if ismember(trialData(i).trial, trialNum)
        idx = find(trialData(i).trial == trialNum);
        kinematicData(i).position = updatedTrialData(idx).kinematics_updated.position;
        kinematicData(i).velocity = updatedTrialData(idx).kinematics_updated.velocity;
        kinematicData(i).acceleration = updatedTrialData(idx).kinematics_updated.acceleration;
        kinematicData(i).speed = updatedTrialData(idx).kinematics_updated.speed;
        kinematicData(i).rotatedPosition = updatedTrialData(idx).kinematics_updated.rotatedPosition;
        kinematicData(i).distanceFromCenter = updatedTrialData(idx).kinematics_updated.distanceFromCenter;
        kinematicData(i).distanceFromEndTarget = updatedTrialData(idx).kinematics_updated.distanceFromEndTarget;
    else
        kinematicData(i).position = NaN;
        kinematicData(i).velocity = NaN;
        kinematicData(i).acceleration = NaN;
        kinematicData(i).speed = NaN;
        kinematicData(i).rotatedPosition = NaN;
        kinematicData(i).distanceFromCenter = NaN;
        kinematicData(i).distanceFromEndTarget = NaN;
    end
end
filename = strrep(kinfile.name, "KIN", "updatedKIN");
save(dataFolder+filename, "kinematicData");

% synchronize the columns in updatedTrialData and behabMetaData
for i = 1:length(trialData)
    behavData(i).trial = trialData(i).trial;
    if ismember(trialData(i).trial, trialNum)
        idx = find(trialData(i).trial == trialNum);
        behavData(i).trialStatusLabel = updatedTrialData(idx).trialStatusLabel;
        behavData(i).t1_trialStartTime = updatedTrialData(idx).t1_trialStartTime;
        behavData(i).t2_targetOnsetTime = updatedTrialData(idx).t2_targetOnsetTime;
        behavData(i).t3_goCueTime = updatedTrialData(idx).t3_goCueTime;
        behavData(i).t4_reactionTime = updatedTrialData(idx).t4_reactionTime;
        behavData(i).t4_reactionTime_ave = updatedTrialData(idx).t4_reactionTime_ave;
        behavData(i).t4b_reactionTime_20mm = updatedTrialData(idx).t4b_reactionTime_20mm;
        behavData(i).t4c_reactionTime_exit = updatedTrialData(idx).t4c_reactionTime_exit;
        behavData(i).t5_peakSpeedTime = updatedTrialData(idx).t5_peakSpeedTime;
        behavData(i).t5_1_homingStartTime = updatedTrialData(idx).t5_1_homingStartTime;
        behavData(i).t6_reachEndTime = updatedTrialData(idx).t6_reachEndTime;
        behavData(i).t6b_reachEndTime_20mm = updatedTrialData(idx).t6b_reachEndTime_20mm;
        behavData(i).t6c_reachEndTime_entry = updatedTrialData(idx).t6c_reachEndTime_entry;
        behavData(i).t6_1_homingEndTime = updatedTrialData(idx).t6_1_homingEndTime;
        behavData(i).t7_postReachStateTime = updatedTrialData(idx).t7_postReachStateTime;
        behavData(i).t8_trialEndTime = updatedTrialData(idx).t8_trialEndTime;
        behavData(i).centerExitSpeed = updatedTrialData(idx).centerExitSpeed;
        behavData(i).extraMetrics.exitMaxSpeed = updatedTrialData(idx).extraMetrics.exitMaxSpeed;
        behavData(i).extraMetrics.exitMaxSpeedOnAxis = updatedTrialData(idx).extraMetrics.exitMaxSpeedOnAxis;
        behavData(i).peakSpeed = updatedTrialData(idx).peakSpeed;
        behavData(i).extraMetrics.avgReachSpeed = updatedTrialData(idx).extraMetrics.avgReachSpeed;
        behavData(i).extraMetrics.timeInTarget = updatedTrialData(idx).extraMetrics.timeInTarget;
        behavData(i).extraMetrics.distInTarget = updatedTrialData(idx).extraMetrics.distInTarget;
        behavData(i).exitTime = updatedTrialData(idx).exitTime;
        behavData(i).midreachTime = updatedTrialData(idx).midreachTime;
        behavData(i).homingTime = updatedTrialData(idx).homingTime;
    else
        % if the trial is not valid, store NaN
        behavData(i).trialStatusLabel = NaN;
        behavData(i).t1_trialStartTime = NaN;
        behavData(i).t2_targetOnsetTime = NaN;
        behavData(i).t3_goCueTime = NaN;
        behavData(i).t4_reactionTime = NaN;
        behavData(i).t4_reactionTime_ave = NaN;
        behavData(i).t4b_reactionTime_20mm = NaN;
        behavData(i).t4c_reactionTime_exit = NaN;
        behavData(i).t5_peakSpeedTime = NaN;
        behavData(i).t5_1_homingStartTime = NaN;
        behavData(i).t6_reachEndTime = NaN;
        behavData(i).t6b_reachEndTime_20mm = NaN;
        behavData(i).t6c_reachEndTime_entry = NaN;
        behavData(i).t6_1_homingEndTime = NaN;
        behavData(i).t7_postReachStateTime = NaN;
        behavData(i).t8_trialEndTime = NaN;
        behavData(i).centerExitSpeed = NaN;
        behavData(i).extraMetrics.exitMaxSpeed = NaN;
        behavData(i).extraMetrics.exitMaxSpeedOnAxis = NaN;
        behavData(i).peakSpeed = NaN;
        behavData(i).extraMetrics.avgReachSpeed = NaN;
        behavData(i).extraMetrics.timeInTarget = NaN;
        behavData(i).extraMetrics.distInTarget = NaN;
        behavData(i).exitTime = NaN;
        behavData(i).midreachTime = NaN;
        behavData(i).homingTime = NaN;
    end
end

% change kinfile.name: replace KIN with BEH
filename = strrep(kinfile.name, "KIN", "BEH");
save(dataFolder+filename, "behavData");

end