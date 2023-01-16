function [trialData] = p5_updateBehavior(trialData,taskInfo)
% This function gets the specific failure mode (or success label) for each
% trial and adds many fields to the kinematics. These include event timings
% (which replaces state transitions), single-value metrics*, ballistic
% reach endpoint predictions, and rotated kinematics. 

% We need the numbers associated with each state
if strcmp(taskInfo.subjectName,'Earl')
    centerAcqState = 2;
    targetOnsetState = 3;
    targetOnPostReactionState = 4; % state after target onset
    goCueState = 5;
    reachState = 6;
    targHoldState = 7;
    rewardState = 8;
    failStates = [9 10];
    failTime = 750;
    targetHoldTime = 400;
    minPS = 0.1; % minimum peak speed for a trial to be considered a reach
elseif strcmp(taskInfo.subjectName,'Prez')
    centerAcqState = 2;
    targetOnsetState = 3;
    targetOnPostReactionState = 3; % state after target onset if needed (or if only one target onset state, just use that)
    goCueState = 4;
    reachState = 4;
    targHoldState = 5;
    rewardState = 7;
    failStates = [8 9];
    failTime = taskInfo.failTime;
    targetHoldTime = taskInfo.targetHoldTime;
    minPS = 0.075;
elseif strcmp(taskInfo.subjectName,'Rocky')
    centerAcqState = 2;
    targetOnsetState = 3;
    targetOnPostReactionState = 3; % state after target onset if needed (or if only one target onset state, just use that)
    goCueState = 4;
    reachState = 5;
    targHoldState = 6;
    rewardState = 7;
    failStates = [10 11 12 13];
    failTime = 667; %taskInfo.failTime
    targetHoldTime = taskInfo.targetHoldTime;
    minPS = 0.1;
end

% Parameters for homing time calc. The behavioral work shows using [1/4 1/3
% 1/2], [0:0.5:2], and [0:50:250] work for these parameters in any combo
% and get basically the same result (respectively); we select this combo to
% best keep Jackpot trials / capture the essence of "homing time"
homingStartDistanceRatio = 1/3; % times target distance - how close to the target is the cursor when homing time begins? e.g. 1/3 = targDist/3 mm away from target entry
homingStopDistanceBuffer = 1; % mm - how far from end target to stop homing time calculation
homingStopTimeBuffer = 150; % ms - how many ms beyond fail time to allow for homing time calc


% Get alll the important timepoints for the reach. the "t#" preceding means
% it's w.r.t. the trial and not a raw quantity (e.g. reactionTime here is
% not true reaction time; it would be RT = t4_reactionTimes-t3_goCueTimes)
t1_trialStartTimes = nan(length(trialData),1);       % what is the first timepoint for this trial?
t2_targetOnsetTimes = nan(length(trialData),1);      % when did the target turn on?
t3_goCueTimes = nan(length(trialData),1);            % when did the go cue happen?
t4_reactionTimes = nan(length(trialData),1);         % when was 20% of peak speed achieved?
t4_reactionTimes_ave = nan(length(trialData),1);     % when was 20% of average peak speed achieved?
t4b_reactionTimes_20mm = nan(length(trialData),1);   % when did speed rise to 0.02m/s = 20 mm/s?
t4c_reactionTimes_exit = nan(length(trialData),1);   % when did the cursor exit the center target?
t5_peakSpeedTimes = nan(length(trialData),1);        % when was peak speed achieved?
t5_1_homingStartTimes = nan(length(trialData),1);    % when does homing time begin?
t6_reachEndTimes = nan(length(trialData),1);         % when was 10% of peak speed achieved (after peak speed)?
t6b_reachEndTimes_20mm = nan(length(trialData),1);   % when did speed fall to 0.02m/s = 20 mm/s (roughly 5% of avg peak speed)
t6c_reachEndTimes_entry = nan(length(trialData),1);  % when did the cursor land in the target?
t6_1_homingEndTimes = nan(length(trialData),1);      % when did homing time end?
t7_postReachStateTimes = nan(length(trialData),1);   % when did the state change from reach to either failure or target hold?
t8_trialEndTimes = nan(length(trialData),1);         % when did reward/failure occur?

bSpeed = 0.02; % speed to use for reaction time and reach end time B

% other important things
kinematics_updated = [];    % kinematics with new additions and aligned and stuff
eyeData_updated = [];       % eye data with new additions and aligned and stuff
centerExitSpeed = nan(length(trialData),1);       % what was the velocity upon exiting the center target?
exitMaxSpeed = nan(length(trialData),1);          % Around exit (+/- 100ms), what was the max speed?
exitMaxSpeedOnAxis = nan(length(trialData),1);    % Around exit (+/- 100ms), what's the max speed towards the target?
peakSpeed = nan(length(trialData),1);             % what was the maximum speed in the reach?
avgReachSpeed = nan(length(trialData),1);         % what was the overall average speed during the reach?
timeInTarget = nan(length(trialData),1);          % how long was the cursor in the end target?
distInTarget = nan(length(trialData),1);          % how much distance did the cursor cover in the target?
trialStatusLabels = nan(length(trialData),1);     % what was the outcome of the trial? All possible outcomes are listed above
goodEyeData = nan(length(trialData),1);           % is the eye data good for this trial?
exitTime = nan(length(trialData),1);              % how long did it take for the animal to exit the target?
midreachTime = nan(length(trialData),1);          % how long did it take the animal to get within 1/3 of the distance to the target?
homingTime = nan(length(trialData),1);            % how long did the animal spend approaching the target?

% For reference
directions = unique([trialData.directionLabel]);
ndirections = length(directions);

for i = 1:length(trialData)
    % get the day (and display text if it changed)
    % set parameters based on the day
    center = trialData(i).centerTarget(1:2);
    startTargRadius = trialData(i).centerTarget(4);
    endTargRadius = trialData(i).reachTarget(4);
    targetDistance = norm(trialData(i).reachTarget(1:2)-trialData(i).centerTarget(1:2));
    curStateTrans = trialData(i).stateTable;
    if ~any(curStateTrans(1,:)==targetOnPostReactionState)
        throw(MException('Er wtf this trial should have been removed? The target was barely presented before failure...'))
    end
    
    % Get overall kinematic stuff we need
    kins = trialData(i).handKinematics;
    time = trialData(i).time;
    posSign = [1 1];
    pos = (kins.position(:,1:2)-center).*posSign; % The negative sign is because the y-axis is flipped (due to the mirror)
    vel = kins.velocity(:,1:2).*posSign;
    acc = kins.acceleration(:,1:2).*posSign;
    targAngle = 45*(trialData(i).directionLabel - 1);
    rotMat = [cosd(targAngle) -sind(targAngle) ; sind(targAngle) cosd(targAngle)];
    rotPos = pos*rotMat;
    distFromCenter = sqrt(pos(:,1).^2+pos(:,2).^2);
    distFromEndTarg = sqrt((targetDistance-rotPos(:,1)).^2+(0-rotPos(:,2)).^2);
    speed = sqrt(sum(vel.^2,2));
    
    % Update eye position
%     if trialData(i).goodEyeData && max(size(trialData(i).eyeData.position)) > 5 % not a short/nan sequence
%         if strcmp(taskInfo.subjectName,'Rocky') % Rocky's signs are flipped already
%             eyeSign = [1 1];
%         else % other animals haven't had eye data flipped yet
%             eyeSign = [1 -1];
%         end
%         eyePos = (trialData(i).eyeData.position(:,1:2)-center).*eyeSign;
%         pupil = trialData(i).eyeData.pupil;
%         goodEyeData(i) = 1;
%     else
%         eyePos = nan(length(time),2);
%         pupil = nan(length(time),1);
%         goodEyeData(i) = 0;
%     end
    
    % Now timings and such. First, let's find when different trial checkpoints happen.
    centerAcqStateInd = find(curStateTrans(1,:)==centerAcqState,1,'last'); % 'last' bc sometimes he quickly enters/exits before target onset, meaning multiple states 1-2
    targetOnsetStateInd = find(curStateTrans(1,:)==targetOnsetState,1);
    goCueStateInd = find(curStateTrans(1,:)==goCueState,1);
    postReachStateInd = find(curStateTrans(1,:)>reachState,1);
    if isempty(centerAcqStateInd) || isempty(targetOnsetStateInd) || isempty(postReachStateInd) % we should have removed all unattempted trials already
        error('WTF we shouldnt have any of these trials anymore')
    else
        t1_trialStartTimes(i) = double(curStateTrans(2,centerAcqStateInd));
        t2_targetOnsetTimes(i) = double(curStateTrans(2,targetOnsetStateInd));
        t7_postReachStateTimes(i) = double(curStateTrans(2,postReachStateInd));
    end
    if isempty(goCueStateInd) % fail before go cue
        t3_goCueTimes(i) = nan;
    else
        t3_goCueTimes(i) = double(curStateTrans(2,goCueStateInd));
    end
    
    
    % Next, we'll get timing/speed associated with exiting the target
    upperExitInd = find((distFromCenter > startTargRadius) & (time>=t2_targetOnsetTimes(i)), 1);
    lowerExitInd = upperExitInd-1;
    if isempty(upperExitInd)
        error('Couldnt find target exit index; this shouldnt happen!')
    else
        t4c_reactionTimes_exit(i) = interp1(distFromCenter(lowerExitInd:upperExitInd),time(lowerExitInd:upperExitInd),startTargRadius);
        exitTime(i) = t4c_reactionTimes_exit(i)-t3_goCueTimes(i);
        centerExitSpeed(i) = interp1(time(lowerExitInd:upperExitInd),speed(lowerExitInd:upperExitInd),t4c_reactionTimes_exit(i));
        buffer = 100/mean(diff(time)); % 100ms
        exitMaxSpeed(i) = max(speed(lowerExitInd-buffer:min(lowerExitInd+buffer,length(speed))));
        exitMaxSpeedOnAxis(i) = max(vel(lowerExitInd-buffer:min(lowerExitInd+buffer,length(speed)),:)*rotMat*[1;0]);
    end
    
    % Next, we'll get peak speed and timings associated with it (shouldn't
    % exist for delay failure trials), including reach end times
    if ~ismember(reachState,curStateTrans(1,:)) % failed before reach
        peakSpeed(i) = nan;
    else
        reachInds = time>=t3_goCueTimes(i) & time<=t7_postReachStateTimes(i);
        if max(speed.*reachInds) < 1
            [peakSpeed(i),peakSpeedInd] = max(speed.*reachInds);
        end
    end
    if ~(~isnan(peakSpeed(i)) && peakSpeed(i) >= minPS) % if it failed before reach or peak speed is too low, stop here
        avgReachSpeed(i) = nan;
    else
        t5_peakSpeedTimes(i) = time(peakSpeedInd);
        lowerReactionInd = find(time<=t5_peakSpeedTimes(i) & speed <= 0.2*peakSpeed(i),1,'last'); % 20% of PS is one way to show RT
        if ~isempty(lowerReactionInd)
            indsToCheck = lowerReactionInd+[0 1];
            if max(speed(indsToCheck)) ~= Inf
                t4_reactionTimes(i) = interp1(speed(indsToCheck),time(indsToCheck),0.2*peakSpeed(i));
            end
        end
        lowerRxn2Ind = find(time<=t5_peakSpeedTimes(i) & speed <= bSpeed,1,'last'); % Using a set value instead of something based on peak speed
        if ~isempty(lowerRxn2Ind) 
            indsToCheck = lowerRxn2Ind+[0 1];
            t4b_reactionTimes_20mm(i) = interp1(speed(indsToCheck),time(indsToCheck),bSpeed);
        end
        upperReachEndInd = find(time>=t5_peakSpeedTimes(i) & speed < 0.1*peakSpeed(i),1);
        if ~isempty(upperReachEndInd)
            indsToCheck = upperReachEndInd-1:upperReachEndInd;
            if max(speed(indsToCheck)) ~= Inf
                t6_reachEndTimes(i) = interp1(speed(indsToCheck),time(indsToCheck),0.1*peakSpeed(i)); % 10% of peak speed is one way to show reach end time
            end
        end
        upperSlowInd = find(time>=t5_peakSpeedTimes(i) & speed < bSpeed,1);
        if ~isempty(upperSlowInd) 
            indsToCheck = upperSlowInd-1:upperSlowInd;
            t6b_reachEndTimes_20mm(i) = interp1(speed(indsToCheck),time(indsToCheck), bSpeed); % Using a set value instead of something based on peak speed
        end
        upperEntryInd = find(time>=t5_peakSpeedTimes(i) & distFromEndTarg <= endTargRadius,1);
        if ~isempty(upperEntryInd)
            indsToCheck = upperEntryInd-1:upperEntryInd;
            t6c_reachEndTimes_entry(i) = interp1(distFromEndTarg(indsToCheck),time(indsToCheck),endTargRadius); % target entry is another way to show reach end time
        end
        avgReachSpeed(i) = mean(speed(time >= t4b_reactionTimes_20mm(i) & time <= t6c_reachEndTimes_entry(i)));
        
        % We also want to calculate midreach and homing times
        homingStartDistance = (targetDistance-startTargRadius-endTargRadius)*homingStartDistanceRatio+endTargRadius;
        homingEndDistance = homingStopDistanceBuffer+endTargRadius;
        maxHomingEndTime = failTime+homingStopTimeBuffer;
        upperHomingStartInd = find((distFromEndTarg<=homingStartDistance) ...
            & (time>=t4c_reactionTimes_exit(i)) ...
            & (time<=t3_goCueTimes(i)+maxHomingEndTime),1);
        if ~isempty(upperHomingStartInd)
            indsToCheck = upperHomingStartInd-[1 0];
            t5_1_homingStartTimes(i) = interp1(distFromEndTarg(indsToCheck),time(indsToCheck),homingStartDistance);
            midreachTime(i) = t5_1_homingStartTimes(i)-t4c_reactionTimes_exit(i);
            upperHomingEndInd = find((distFromEndTarg<=homingEndDistance) ...
                & (time>=t5_1_homingStartTimes(i)) ...
                & (time<=t3_goCueTimes(i)+maxHomingEndTime),1);
            if ~isempty(upperHomingEndInd)
                indsToCheck = upperHomingEndInd-[1 0];
                t6_1_homingEndTimes(i) = interp1(distFromEndTarg(indsToCheck),time(indsToCheck),homingEndDistance);
                homingTime(i) = t6_1_homingEndTimes(i)-t5_1_homingStartTimes(i);
            end
        end
    end
    
    % And finally, get trial "end time", which we'll note as either time
    % time of reward or failure
    t8_trialEndTimes(i) = curStateTrans(2,find(ismember(curStateTrans(1,:),[rewardState failStates]),1,'first'));
    
    % If the cursor reached the target in time, we can calculate how
    % long and how much distance was covered while in the target
    if ismember(targHoldState,curStateTrans(1,:))
        timeInTarget(i) = t8_trialEndTimes(i)-t6c_reachEndTimes_entry(i);
        indsInTarget = time>=t6c_reachEndTimes_entry(i) & time<=t8_trialEndTimes(i);
        distInTarget(i) = sum(sqrt(sum(diff(pos(indsInTarget,:)).^2,2)));
    end
    
    
    % We calculate a predicted ballistic reach endpoint. We can do this by
    % mirroring the velocity traces up to peak speed and finding the
    % endpoint. An equal way to do this is to double the displacement from
    % start of reach to peak speed; this is easier, so we'll do that.
    if isnan(peakSpeed(i)) || peakSpeed(i) < minPS || peakSpeed(i) > 1 % if a no reach occurred, nan it (or if it's a crazy trial)
        ballisticEndpointPrediction = nan(1,2);
        rotBallisticEndpointPrediction = nan(1,2);
    else
        firstInd = round((t4b_reactionTimes_20mm(i)-50)/mean(diff(time))); % we use 50ms before the reach begins as our starting point
        dispToPS = diff(pos([firstInd peakSpeedInd],:));
        ballisticEndpointPrediction = dispToPS*2+pos(firstInd,:);
        rotBallisticEndpointPrediction = ballisticEndpointPrediction*rotMat;
    end
    
    
    % Now we identify the status of the trial as just success or failure;
    % we need the statistics for successful trials for some failure
    % identification, so for now we'll just label success (1) or fail (0)
    trialStatusLabels(i) = ismember(rewardState,curStateTrans(1,:));
    
    % Now we just have to assign all of the stuff for a given trial to it's
    % appropriate substructure (if it has one)...
    kinematics_updated(i).position = pos;
    kinematics_updated(i).velocity = vel;
    kinematics_updated(i).acceleration = acc;
    kinematics_updated(i).speed = speed;
    kinematics_updated(i).rotatedPosition = rotPos;
    kinematics_updated(i).distanceFromCenter = distFromCenter;
    kinematics_updated(i).distanceFromEndTarget = distFromEndTarg;
    kinematics_updated(i).ballisticEndpointPrediction = ballisticEndpointPrediction;
    kinematics_updated(i).rotatedBallisticEndpointPrediction = rotBallisticEndpointPrediction;
    
%     eyeData_updated(i).position = eyePos;
%     eyeData_updated(i).pupil = pupil;
    
%     if strcmp(taskInfo.subjectName,'Rocky') % Rocky has EMG and heart for some extra processing
%         heartData_updated(i) = trialData(i).heartData_EKG;
%         goodHeartData(i) = ~any(isnan(heartData_updated(i).heartRate));
%         EMGData_updated(i) = trialData(i).EMGData;
%         goodEMGData(i,:) = [1 1 1 1 1];
%         % ADAM HERE we need to check EMG for drops again here or something.
%         % Update the "good Emg data" field
%     end
end; clear i
% disp(['Completed behavior updating for day ' num2str(day)])

for i = 1:length(trialData)
    if ~isnan(t4_reactionTimes(i))
        kins = trialData(i).handKinematics;
        time = trialData(i).time;
        vel = kins.velocity(:,1:2);
        speed = sqrt(sum(vel.^2,2));
        lowerReactionInd = find(time<=t5_peakSpeedTimes(i) & speed <= 0.2*mean(peakSpeed(~isnan(peakSpeed))),1,'last'); % 20% of PS is one way to show RT
        if ~isempty(lowerReactionInd)
            indsToCheck = lowerReactionInd+[0 1];
            t4_reactionTimes_ave(i) = interp1(speed(indsToCheck),time(indsToCheck),0.2*mean(peakSpeed(~isnan(peakSpeed))));
        end
    end
end



% Now, we go through the failed trials and give them a failure label:
%   0 = Quitout: an active reach away from the direction of the reach target
% -11 = False Start: a reach attempt was made towards the targetbefore the go cue had been presented
% -12 = Delay Drift: hand lazily drifted out of the center, sometimes towards the target ("cheat") others not
% -20 = Wild: reach is too far off the course to call it either an overshoot or undershoot
% -21 = No Attempt: no attempt (or a very lazy one) was made, indicated by a very low or super late peak speed
% -22 = Overshoot (type 1: miss): reach missed target
% -23 = Undershoot (type 1: inaccurate): reach movement ended before target
% -24 = Undershoot (type 2: slow): reach was still ongoing when time expired, which occurred short of the target
% -31 = Overshoot (type 2: scuff): cursor barely scrapes the target and goes through it
% -32 = Overshoot (type 3: blow-through): reach goes right through the target
% -33 = Early Return: reaching back to center before time is up
% -34 = Target Hold Drift: hand drifts out of the end target
% -35 = Jitter: land at edge of target -> small jitter pulls cursor out
%
% -1x = delay epoch, -2x = reach epoch, -3x = target hold epoch
thresholdMaxExitSpeed = prctile(exitMaxSpeedOnAxis(trialStatusLabels==1),1); % center exit speed above this = false start or quitout
for i = 1:length(trialData)
    if trialStatusLabels(i)~=1
        time = trialData(i).time;
        endTargRadius = trialData(i).reachTarget(4);
        targetDistance = norm(trialData(i).reachTarget(1:2)-trialData(i).centerTarget(1:2));
        curStateTrans = trialData(i).stateTable;
        if ~ismember(reachState,curStateTrans(1,:)) % some form of delay failure
            if exitMaxSpeedOnAxis(i) > thresholdMaxExitSpeed % reached towards target with some speed
                % False Start: a reach attempt was made towards the target
                % before the go cue had been presented
                trialStatusLabels(i) = -11;
            elseif exitMaxSpeed(i) > thresholdMaxExitSpeed % reached away from target with some speed
                % Quitout: an active reach away from the direction of the
                % reach target
                trialStatusLabels(i) = 0;
            else
                % Delay Drift: hand lazily drifted out of the center,
                % sometimes towards the target ("cheat") others not
                trialStatusLabels(i) = -12;
            end
        elseif ~ismember(targHoldState,curStateTrans(1,:)) % target not acquired; reach epoch failure
            goCueInd = find(time>=t3_goCueTimes(i),1);
            failInd = find(time>=t7_postReachStateTimes(i),1);
            rotPos = kinematics_updated(i).rotatedPosition;
            angleToFailRotPos = atan2d(rotPos(failInd,2),rotPos(failInd,1));
            if (peakSpeed(i) < minPS) || ((exitTime(i) < 400) && (t7_postReachStateTimes(i)-t5_peakSpeedTimes(i) < 50)) % either super low peak speed or very late peak speed
                % No Attempt: no attempt (or a very lazy one) was made,
                % indicated by a very low or super late peak speed
                trialStatusLabels(i) = -21;
            elseif abs(angleToFailRotPos) > 60 % cursor @failure didn't end up within 60deg of vector from center to target
                % Quitout: an active reach away from the direction of the
                % reach target
                trialStatusLabels(i) = 0;
            elseif peakSpeed(i) > 1 % some odd cursor jump that didn't get caught; remove
                trialStatusLabels(i) = 0;
            elseif max(abs(rotPos(goCueInd:failInd,2))) > 25 % reach had component way off of the line to the target
                % Wild: reach is too far off the course to call it either
                % an overshoot or undershoot
                trialStatusLabels(i) = -20;
            elseif max(rotPos(goCueInd:failInd,1)) > targetDistance % reached past center of end target
                % Overshoot (type 1: miss): reach missed target
                trialStatusLabels(i) = -22;
            elseif t6b_reachEndTimes_20mm(i) < t7_postReachStateTimes(i) % reach ended (came to very low speed) outside of target
                % Undershoot (type 1: inaccurate): reach movement ended before target
                trialStatusLabels(i) = -23;
            else % mid-reach when time expired
                % Undershoot (type 2: slow): reach was still ongoing when time
                % expired, which occurred short of the target
                trialStatusLabels(i) = -24;
            end
        else % target acquired but left; target hold epoch failure
            enterInd = find(time>=t7_postReachStateTimes(i),1);
            exitInd = find(time>=t8_trialEndTimes(i),1);
            rotPos = kinematics_updated(i).rotatedPosition;
            reachTargExitSpeed = kinematics_updated(i).speed(exitInd);
            angleOfTargExitRotPos = atan2d(rotPos(exitInd,2),rotPos(exitInd,1)-targetDistance);
            if exitInd-enterInd <= 1 % weird bug; for some reason the trial is a failure immediately on entry, even if it doesn't look like a scuff...throw it out
                error(['Bad trial, i = ' num2str(i) ' (immediate entry failure) - remove this and re-run'])
            elseif (timeInTarget(i) < 100) && (distInTarget(i) < endTargRadius) % not long spent in target, and not much distance covered in it
                % Overshoot (type 2: scuff): cursor barely scrapes the
                % target and goes through it
                trialStatusLabels(i) = -31;
            elseif (timeInTarget(i) < targetHoldTime/2) && (distInTarget(i) > endTargRadius)
                % Overshoot (type 3: blow-through): reach goes right
                % through the target
                trialStatusLabels(i) = -32;
            elseif (reachTargExitSpeed > 0.1) && (abs(angleOfTargExitRotPos) > 90) % somewhat fast exit, somewhat towards center
                % Early Return: reaching back to center before time is up
                trialStatusLabels(i) = -33;
            elseif distInTarget(i) > endTargRadius % doesn't blow through target, but drifts out
                % Target Hold Drift: hand drifts out of the end target
                trialStatusLabels(i) = -34;
            else % not much distance in target
                % Jitter: land at edge of target -> small jitter pulls
                % cursor out
                trialStatusLabels(i) = -35;
            end
        end
    end
end; clear i
disp('Completed failure mode evaluation')


% ...and finally, assign these all to trialData! We'll put stuff that isn't
% as important in an "extraMetrics" structure
trialData = rmfield(trialData,'handKinematics');
for i = 1:length(trialData)
    trialData(i).kinematics_updated = kinematics_updated(i);
%     trialData(i).eyeData_updated = eyeData_updated(i);
%     trialData(i).goodEyeData = goodEyeData(i);
%     if strcmp(taskInfo.subjectName,'Rocky')
%         trialData(i).heartData_updated = heartData_updated(i);
%         trialData(i).goodHeartData = goodHeartData(i);
%         trialData(i).EMGData_updated = EMGData_updated(i);
%         trialData(i).goodEMGData = goodEMGData(i,:);
%     end
    trialData(i).trialStatusLabel = trialStatusLabels(i);
    trialData(i).t1_trialStartTime = t1_trialStartTimes(i);
    trialData(i).t2_targetOnsetTime = t2_targetOnsetTimes(i);
    trialData(i).t3_goCueTime = t3_goCueTimes(i);
    trialData(i).t4_reactionTime = t4_reactionTimes(i);
    trialData(i).t4_reactionTime_ave = t4_reactionTimes_ave(i);
    trialData(i).t4b_reactionTime_20mm = t4b_reactionTimes_20mm(i);
    trialData(i).t4c_reactionTime_exit = t4c_reactionTimes_exit(i);
    trialData(i).t5_peakSpeedTime = t5_peakSpeedTimes(i);
    trialData(i).t5_1_homingStartTime = t5_1_homingStartTimes(i);
    trialData(i).t6_reachEndTime = t6_reachEndTimes(i);
    trialData(i).t6b_reachEndTime_20mm = t6b_reachEndTimes_20mm(i);
    trialData(i).t6c_reachEndTime_entry = t6c_reachEndTimes_entry(i);
    trialData(i).t6_1_homingEndTime = t6_1_homingEndTimes(i);
    trialData(i).t7_postReachStateTime = t7_postReachStateTimes(i);
    trialData(i).t8_trialEndTime = t8_trialEndTimes(i);
    trialData(i).centerExitSpeed = centerExitSpeed(i);
    trialData(i).extraMetrics.exitMaxSpeed = exitMaxSpeed(i);
    trialData(i).extraMetrics.exitMaxSpeedOnAxis = exitMaxSpeedOnAxis(i);
    trialData(i).peakSpeed = peakSpeed(i);
    trialData(i).extraMetrics.avgReachSpeed = avgReachSpeed(i);
    trialData(i).extraMetrics.timeInTarget = timeInTarget(i);
    trialData(i).extraMetrics.distInTarget = distInTarget(i);
    trialData(i).exitTime = exitTime(i);
    trialData(i).midreachTime = midreachTime(i);
    trialData(i).homingTime = homingTime(i);
end; clear i

% % If desired, order the fields
% fieldNameOrder = sort(fieldnames(trialData(1)));
% if strcmp(taskInfo.subjectName,'Rocky')
%     fieldNameOrder = {'day','trial','directionLabel','rewardLabel','trialStatusLabel','goodEyeData','goodHeartData','goodEMGData',... % labels
%     'photoTransitions',... % state transition timings
%     'time','kinematics_updated','eyeData_updated','heartData_updated','EMGData_updated','neuralData',... % time series data
%     't1_trialStartTime','t2_targetOnsetTime','t3_goCueTime','t4_reactionTime',... % timings
%     't4b_reactionTime_20mm','t4c_reactionTime_exit','t5_peakSpeedTime','t5_1_homingStartTime',...
%     't6_reachEndTime','t6b_reachEndTime_20mm','t6c_reachEndTime_entry','t6_1_homingEndTime',...
%     't7_postReachStateTime','t8_trialEndTime','exitTime','midreachTime','homingTime',... % also metrics for timings
%     'centerExitSpeed','peakSpeed','extraMetrics',... % everything else
%     };
% else
%     fieldNameOrder = {'day','trial','directionLabel','rewardLabel','trialStatusLabel','goodEyeData',... % labels
%     'photoTransitions',... % state transition timings
%     'time','kinematics_updated','eyeData_updated','neuralData',... % time series data
%     't1_trialStartTime','t2_targetOnsetTime','t3_goCueTime','t4_reactionTime',... % timings
%     't4b_reactionTime_20mm','t4c_reactionTime_exit','t5_peakSpeedTime','t5_1_homingStartTime',...
%     't6_reachEndTime','t6b_reachEndTime_20mm','t6c_reachEndTime_entry','t6_1_homingEndTime',...
%     't7_postReachStateTime','t8_trialEndTime','exitTime','midreachTime','homingTime',... % also metrics for timings
%     'centerExitSpeed','peakSpeed','extraMetrics',... % everything else
%     };
% end
% trialData = orderfields(trialData, fieldNameOrder);

end

