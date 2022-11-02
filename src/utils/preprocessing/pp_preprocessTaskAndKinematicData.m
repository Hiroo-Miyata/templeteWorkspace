function [trialData,taskInfo,kinematicData] = pp_preprocessTaskAndKinematicData(taskDataFolders,rewardNames,synchInfo)
% This fucntion preprocesses the behavioral data for each trial of Rocky's
% delayed center out task for further use. Specifically, it does 3
% things:
% (1) Preprocesses position data to get velocity/acceleration, smooths
% them, and interpolates to 1000Hz. It will also convert everything into
% mm centered on the center target
% (2) Extracts trial labels for reach target direction, reward size, and
% whether or not it's a catch trial
% (3) Converts everything into task coordinates - mm, and centered on the
% task center
%
% The output is the following structures:
% - trialData: structure with the following fields:
%     - taskSynchTrialTime: scalar, time of trial synch pulse (center target appearance) - this is what you'd use to align to the analog data I believe.
%     - trial: scalar, trial #
%     - attempted: 0/1, whether or not he attempted the trial (should all be 1)
%     - trialStatus: 0/1, 0 = failure, 1 = success
%     - stateTable: [2 x nstates] array of when the state transitions in this trial happened. Top row is state number, second row is the time (in ms!) from the         start of the trial when the state occurred
%         it has been set such that the center target appearance (state 1) is time 0
%     - centerTarget: [X Y Z r] of the center target in mm (r is radius)
%     - reachTarget: [X Y Z r] of the reach target in mm
%     - directionLabel: integer on [1:8] indicating reach target direction, where 1 = 0deg, 2 = 45deg, 3 = 90deg, 4 = 135deg, 5 = 180deg, 6 = 225deg, 7 = 270deg, 8 = 315deg
%     - rewardLabel: integer on [1:4] indicating reward cue for the trial. 1 = Small (0uL), 2 = Medium (210uL), 3 = Large (420uL), 4 = Jackpot (2100uL)
%     - prevStatusLabel: label of success/fail for previous trial (0 = fail, 1 = success). 1st trial of the day is just set to 0.
%     - prevDirectionLabel: previous direction label; set to 1 for first trial each day
%     - prevRewardLabel: previous reward label; set to 1 (I think?) for first trial each day
%     - time: [ntime x 1] relative time for this trial, where 0 = center target appearance
% - kinematicData: structure with the following fields
%     - position/velocity/acceleration: [ntime x 3] matrices with the kinematics interpolated to 1KHz for this trial in [X Y Z]
% - taskInfo: structure with information about task parameters. Contains many fields, but some important ones include:
%     - centerHoldTime: Amount of time in ms the center is held (state 2) before the target shows up (state 3)
%     - targetHoldTime: Amount of time in ms that the reach target must be held (state 6) until the reward state (state 7) is entered
%     - lim_reachTime: scalar, maximum time allowed after go cue (state 4) until the trial is deemed a failure
%     - stateInfo: name and # associated with each state in stateTable
%     - ppInfo_KIN: contains parameters used for kinematic processing
%
% Assumes the Target.m class is somewhere in the path
%
% Adam Smoulder, 11/2/22
%
% NOTE: The top section here will need edited for datasets with more than
% one task data folder (i.e. expts with task restarts between blocks). To
% do this, you just need to loop over the taskDataFolders when assembling
% instead of just accessing element {1}
% trialData. 

%% First, load all of the trial data
disp('Beginning preprocessing for task and kinematic data')
disp(taskDataFolders)
allFiles = dir(taskDataFolders{1});
nfiles = size(allFiles,1);
trialCount = 0;
goodInds = [];
for i = 1:nfiles
    curFileName = allFiles(i).name;
    if length(curFileName)>=14 && (strcmp(curFileName(1:5),'Trial') && strcmp(curFileName(end-3:end),'.mat')) % If a trial file
        trialCount = trialCount+1;
        curTrialData = load([taskDataFolders{1} curFileName]);
        if trialCount==1
            trialData = curTrialData.trialData;
        else
            trialData(end+1) = curTrialData.trialData;
        end
        if mod(trialCount,100)==0, disp(['Loaded ' num2str(trialCount) 'th trial']); end
        
        % We'll also identify the config file indices for later
    elseif length(curFileName)>=14 && strcmp(curFileName(1:19),'TrialDataConfigFile')
        trialFileInd = i;
    elseif length(curFileName)>=14 &&strcmp(curFileName(1:21),'ReachTargetConfigFile')
        targetFileInd = i;
    end
end; clear i
trialData = trialData';
ntrials = length(trialData);

%  Load taskInfo and add parameters if they're missing; even if multiple 
%  sessions are combined, for ease, we're going to assume taskInfo from the
%  last session is relevant
load([taskDataFolders{1} 'taskInfo.mat'])
displayParams = taskInfo.displayParams;
if ~any(contains(fieldnames(displayParams),'AFFECTOR_CENTER'))
    disp('Adding default AFFECTOR_CENTER')
    taskInfo.displayParams.AFFECTOR_CENTER = [0.0845 0.0084 -1.5];
end
if ~any(contains(fieldnames(displayParams),'AFF_TO_SCREEN_SF'))
    disp('Adding default AFF_TO_SCREEN_SF')
    taskInfo.displayParams.AFF_TO_SCREEN_SF = 1000;
end
if ~any(contains(fieldnames(displayParams),'COMFORT_SF'))
    disp('Adding default COMFORT_SF')
    taskInfo.displayParams.COMFORT_SF = 1;
end


%% Go through each trial and extract/filter all things needed

% Set parameters
targetAngles = 0:45:315; disp('ASSUMING DEFAULT 8 TARGET LOCATIONS')
initialFs = 60; % nominal frequency of the hand data
finalFs = 1000; % final desired frequency to interpolate to
fc = 10;        % cutoff frequency of LPF for kinematics
order = 4/2;    % order (/2) of LPF for kinematics; filtfilt means the order is applied twice, hence /2
[bLPF,aLPF] = butter(order,fc/(initialFs/2),'low'); % make the filter
trialStartState = 1;

% In case it's the choice task, we need a few additional things
target1RewardState = 7;
target2RewardState = 9;
reachTargetNames = {taskInfo.reachTargets.Name}';


% Run it!
nunattempted = 0; % number of unattempted trials
ncursorjump = 0; % number of cursor jump trials
for i = 1:length(trialData)
    if isfield(synchInfo,'badSynchInds') && ismember(i,synchInfo.badSynchInds)
        newTrialData(i).taskSynchTrialTime = nan;
    else
        newTrialData(i).taskSynchTrialTime = synchInfo.taskSynchTrialTimes(i);
    end
    newTrialData(i).trial = trialData(i).trial;
    newTrialData(i).attempted = trialData(i).attempted;
    newTrialData(i).trialStatus = trialData(i).trialStatus;
    newTrialData(i).stateTable = trialData(i).stateTable;
    
    % Depending on if it's choice or main task, we have diff things to get 
    if isfield(trialData(i),'reachTarget') % delayed cout task
        % Get the reach and center target info, properly aligning the y-axis
        newTrialData(i).centerTarget = [trialData(i).centerTarget.Location.*[1 -1]  0  trialData(i).centerTarget.Radius]; % this is the format Batista lab uses
        newTrialData(i).reachTarget = [trialData(i).reachTarget.Location.*[1 -1]  0  trialData(i).reachTarget.Radius]; % this is the format Batista lab uses
        
        % Get the target label from the direction of the target w.r.t. center
        curTargetXY = newTrialData(i).reachTarget(1:2)-newTrialData(i).centerTarget(1:2);
        curTargetAngle = atan2d(curTargetXY(2),curTargetXY(1));
        [~,newTrialData(i).directionLabel] = min(abs(targetAngles-mod(curTargetAngle+360,360)));  % The mod and + is to make it on [0 359] degrees
        
        % Get the reward and catch label from the target name
        reachTargetName = trialData(i).reachTarget.Name;
        newTrialData(i).rewardLabel = find(cellfun(@(x) contains(reachTargetName,x), rewardNames));
        newTrialData(i).catchLabel = contains(reachTargetName,'Catch');
        
    elseif isfield(trialData(i),'reachTarget1') % choice task (has reachTarget2 as well)
        % Get the reach and center target info, properly aligning the y-axis
        newTrialData(i).centerTarget = [trialData(i).centerTarget.Location.*[1 -1]  0  trialData(i).centerTarget.Radius]; % this is the format Batista lab uses
        newTrialData(i).reachTarget1 = [trialData(i).reachTarget1.Location.*[1 -1]  0  trialData(i).reachTarget1.Radius]; % this is the format Batista lab uses
        newTrialData(i).reachTarget2 = [trialData(i).reachTarget2.Location.*[1 -1]  0  trialData(i).reachTarget2.Radius]; % this is the format Batista lab uses
        
        % Get the target label from the direction of the target w.r.t. center
        curTarget1XY = newTrialData(i).reachTarget1(1:2)-newTrialData(i).centerTarget(1:2);
        curTarget1Angle = atan2d(curTarget1XY(2),curTarget1XY(1));
        [~,newTrialData(i).direction1Label] = min(abs(targetAngles-mod(curTarget1Angle+360,360)));  % The mod and + is to make it on [0 359] degrees
        
        % Get the reward and catch label from the target name
        reachTarget1Name = trialData(i).reachTarget1.Name;
        newTrialData(i).reward1Label = find(cellfun(@(x) contains(reachTarget1Name,x), rewardNames));
        newTrialData(i).catch1Label = contains(reachTarget1Name,'Catch');
        
        % Do these both for target 2 as well
        curTarget2XY = newTrialData(i).reachTarget2(1:2)-newTrialData(i).centerTarget(1:2);
        curTarget2Angle = atan2d(curTarget2XY(2),curTarget2XY(1));
        [~,newTrialData(i).direction2Label] = min(abs(targetAngles-mod(curTarget2Angle+360,360))); 
        reachTarget2Name = trialData(i).reachTarget2.Name;
        newTrialData(i).reward2Label = find(cellfun(@(x) contains(reachTarget2Name,x), rewardNames));
        newTrialData(i).catch2Label = contains(reachTarget2Name,'Catch');
        
        % Evaluate which target was chosen
        if any(trialData(i).stateTable(1,:)==target1RewardState)
            newTrialData(i).choiceLabel = 1;
            newTrialData(i).directionLabel = newTrialData(i).direction1Label; 
            newTrialData(i).rewardLabel = newTrialData(i).reward1Label;
            newTrialData(i).catchLabel = newTrialData(i).catch1Label;
        elseif any(trialData(i).stateTable(1,:)==target2RewardState)
            newTrialData(i).choiceLabel = 2;
            newTrialData(i).directionLabel = newTrialData(i).direction2Label;
            newTrialData(i).rewardLabel = newTrialData(i).reward2Label;
            newTrialData(i).catchLabel = newTrialData(i).catch2Label;
        else % no reward; figure it out later
            newTrialData(i).choiceLabel = 0;
            newTrialData(i).directionLabel = 0;
            newTrialData(i).rewardLabel = 0;
            newTrialData(i).catchLabel = 0;
        end
        
    else
        error('Could not find reachTarget or reachTarget1 property on data')
    end
    
    % Get previous parameters too; if previous trial was unattempted, set all to 0
    if i > 1 && newTrialData(i-1).attempted
        newTrialData(i).prevStatusLabel = newTrialData(i-1).trialStatus;
        newTrialData(i).prevDirectionLabel = newTrialData(i-1).directionLabel;
        newTrialData(i).prevRewardLabel = newTrialData(i-1).rewardLabel;
    else
        newTrialData(i).prevStatusLabel = 0;
        newTrialData(i).prevDirectionLabel = 0;
        newTrialData(i).prevRewardLabel = 0;
    end
    
    
    
    % Get the kinematics from the affector and clean them up
    % First, just get the current position and time
    curKinTime = trialData(i).time*1000; % w.r.t. start of trial; in ms
    curPos = (trialData(i).affectorPosition-taskInfo.displayParams.AFFECTOR_CENTER)*taskInfo.displayParams.AFF_TO_SCREEN_SF.*[1 -1 1]; 
    if i == 1
        z_center = median(curPos(:,3)); % the affector center for z is meaningless since it's unused, so set an arbitrary mean
    end
    curPos(:,3) = curPos(:,3)-z_center;
    

    % If it's a failure within < 100ms of target appear, call unattempted
    stateTable = trialData(i).stateTable;
    targetOnsetState = 3;
    targetOnsetTime = stateTable(2,stateTable(1,:)==targetOnsetState);
    if ~isempty(targetOnsetTime) % target acquired
        delayLength = stateTable(2,find(stateTable(1,:)==targetOnsetState)+1)-targetOnsetTime; % time of next state (GC, delay fail, or catch success)
        if delayLength < 100 % <100ms = not enough time to perceive target
            newTrialData(i).attempted = false;
        end
    else % unattempted
        newTrialData(i).attempted = false;
    end
    
    
    % Evaluate if the trial has any cursor jumps
    if newTrialData(i).attempted
        centerAcqState = 2;
        centerAcqTime = stateTable(2,find(stateTable(1,:)==centerAcqState,1,'last')); % technically center acq can go in and out
        validInds = curKinTime >= centerAcqTime;
        validSpeed = sqrt(sum((diff(curPos(validInds,:))./diff(curKinTime(validInds))).^2,2)); % mm/ms = m/s
%         figure; plot(validSpeed,'o-'); set(gcf,'position',[-3839 -222 1920 963]); title([num2str(i), ' max = ' num2str(max(validSpeed))]);
%         pause; close all % use this for debugging
        if any(validSpeed > 1) % far as I've seen, speed is < 0.7-0.8m/s in reality, so > 1 = cursor drop
            newTrialData(i).cursorJump = true;
            ncursorjump = ncursorjump+1;
        else
            newTrialData(i).cursorJump = false;
        end
    else % unattempted; add to counter
        nunattempted = nunattempted+1;
        newTrialData(i).cursorJump = false;
    end
    
    % Our sampling is sometimes uneven, which we can't filter on; but, the
    % data are pretty smooth as is, so we can confidently interpolate to a
    % valid sampling rate. We'll use 60Hz, the nominal resolution
    dt = 1000/initialFs;
    posTime = (curKinTime(1):dt:curKinTime(end))';
    curPos_rs = interp1(curKinTime,curPos,posTime,'spline');
    
    % Filter the data and calculate velocity and acceleration
    curPos_filt = filtfilt(bLPF,aLPF,double(curPos_rs)); % in mm
    curVel_filt = diff(curPos_filt)./diff(posTime); % mm/ms = m/s
    velTime = (posTime(1:end-1)+posTime(2:end))/2; % this velocity approximates that between time indices
    curAcc_filt = diff(curVel_filt)./(diff(velTime)/1000); % (m/s)/(ms/1000) = (m/s)/s = m/s^2
    accTime = posTime(2:end-1); % same as doing the averages for the velocity time
    
    % Interpolate to the new fs and align to center target appearance
    time = (ceil(min(posTime)):1000/finalFs:floor(max(posTime)))';
    kinematicData(i).position = interp1(posTime,curPos_filt,time,'spline');
    kinematicData(i).velocity = interp1(velTime,curVel_filt,time,'spline');
    kinematicData(i).acceleration = interp1(accTime,curAcc_filt,time,'spline');
    centerTargAppearTime = trialData(i).stateTable(2,find(trialData(i).stateTable(1,:)==trialStartState,1));
    newTrialData(i).time = time-centerTargAppearTime;
    newTrialData(i).stateTable(2,:) = newTrialData(i).stateTable(2,:)-centerTargAppearTime;
    
    if mod(i,100)==0, disp(['Processed trial ' num2str(i)]); end
end; clear i
trialData = newTrialData'; clear newTrialData;
disp(['Identified ' num2str(ncursorjump) ' cursor jump trials and ' num2str(nunattempted) ' unattempted trials to be removed']) 

% Get across-trial info too
ppInfo_KIN.originalFs = initialFs;
ppInfo_KIN.newFs = finalFs;
ppInfo_KIN.LPFOrder = order;
ppInfo_KIN.fc = fc;

% clean up taskInfo a bit
taskParams = rmfield(taskInfo,{'subjectName','datetime'});
taskInfo = rmfield(taskInfo,fieldnames(taskParams));
taskInfo.taskParams = taskParams;
taskInfo.ppInfo_KIN = ppInfo_KIN;

end

