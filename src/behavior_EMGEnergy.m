close all; clear all;
dates = ["0407", "0408", "0412", "0413", "0414", "0419", "0420"]; % "0407", "0408", "0412", "0413", "0414", "0415", "0419", "0420"
fileName = "all";
resultDir = "../results/202306w3-results-summary-"+ fileName;
outputFolder = resultDir + "/behavior_EMG/";

%% conbine across days
trialDataAll = struct.empty;
behavDataAll = struct.empty;
EMGDataAll = struct.empty;
EMGMetricsAll = cell(length(dates), 1);
trialNumBegin = 0;
for d = 1:length(dates)
    date = dates(d);
    dataFolder = "../data/preprocessed/2022"+date+"/";
    taskfile = dir(dataFolder+ "*" +date+ "*_TSK_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    taskfile = dir(dataFolder+ "*" +date+ "*_BEH_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    taskfile = dir(dataFolder+ "*" +date+ "*_EMG_*.mat");
    load(taskfile.folder +"/"+ taskfile.name); 

    trialNum = [trialData.trial];
    curInds = trialNum > 24;
    trialData = trialData(curInds);
    behavData = behavData(curInds);
    EMGData = EMGData(curInds);

    for i=1:length(trialData)
        trialData(i).dayLabel = d;
        trialData(i).newTrial = trialData(i).trial + trialNumBegin;
    end

    trialDataAll = cat(1, trialDataAll, trialData); clear trialData
    behavDataAll = cat(2, behavDataAll, behavData); clear behavData
    EMGDataAll = cat(2, EMGDataAll, EMGData); clear EMGData
    EMGMetricsAll{d} = EMGMetrics; clear EMGMetrics
    trialNumBegin = trialNumBegin + max(trialNum);
end; clear taskfile date dataFolder
trialData = trialDataAll; clear trialDataAll
behavData = behavDataAll; clear behavDataAll
EMGData = EMGDataAll; clear EMGDataAll
EMGMetrics = EMGMetricsAll; clear EMGMetricsAll

nmuscles = 5;

%% Use only success trials (effort = reaching + holding)
curInds = false(length(trialData), 1);
for i=(1:length(trialData))
    stateTransition = trialData(i).stateTable;
    if all(ismember([3 4 5 6 7], stateTransition(1,:))) == 1
        curInds(i) = true;
    end
end

trialData = trialData(curInds);
behavData = behavData(curInds);
EMGData = EMGData(curInds);

%% get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels = cat(1, trialData(:).reachTarget); difficultyLabels = difficultyLabels(:, 4); difficultyLabels = (difficultyLabels==9).';
difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
trialStatusLabels = [behavData.trialStatusLabel];
delayTimes = [behavData.delayTime]; delayTimes = round(delayTimes / 50) * 50;
reactionTimes = [behavData.reactionTime];
peakSpeeds = [behavData.peakSpeed];
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
DiffStyle = ["-", ":"];
DelayTimes = 300:50:1050; nDelayTimes = length(DelayTimes);
muscleLabel = EMGMetrics{1}.muscleNames;


energy = zeros(length(trialData),5);
for i=1:length(trialData)
    stateTransition = trialData(i).stateTable;
    %% Calculate EMG TrajMov(1/𝑁𝑚𝑢𝑠𝑐𝑙𝑒𝑠 ∑_(𝑖∈𝑚𝑢𝑠𝑐𝑙𝑒𝑠)▒∫1▒〖‖〖𝐸𝑀𝐺〗_𝑖 ‖^2  𝑑𝑡〗) during movement(from GC to Success)
    GoCueTime = stateTransition(2, find(stateTransition(1, :)==4));
    GoCueTime = find(trialData(i).time == GoCueTime);
    EndTime = stateTransition(2, find(stateTransition(1, :)==7));
    EndTime = find(trialData(i).time == EndTime-75);

    EMGofMovement = movmean(EMGData(i).signal(GoCueTime:EndTime, :), 200); % from GC to Success
    for ch = 1:nmuscles
        meanE = EMGMetrics{1}.normalizedParams(1, ch);
        stdE = EMGMetrics{1}.normalizedParams(2, ch);
%             eachE(ch) = meanE.^2 *size(EMGofMovement, 1) + sumsqr(EMGofMovement(:,ch)) * stdE.^2 ...
%                         + 2*stdE*meanE*sum(EMGofMovement(:,ch));
        energy(i, ch) = meanE *size(EMGofMovement, 1) + stdE*sum(EMGofMovement(:,ch));
        energy(i, ch) = energy(i, ch) * 0.001;
    end
end


%% Plotting
for ch=1:nmuscles
    Y = energy(:, ch)';
    VS_reward_size(Y, rewardLabels, difficultyLabels, "EMG energy: "+muscleLabel(ch)+" (μV)", outputFolder+muscleLabel(ch))
    VS_delay_time(Y, rewardLabels, difficultyLabels, delayTimes, "EMG energy: "+muscleLabel(ch)+" (μV)", outputFolder+muscleLabel(ch))
    VS_direction(Y, rewardLabels, difficultyLabels, directionLabels, "EMG energy: "+muscleLabel(ch)+" (μV)", outputFolder+muscleLabel(ch))
    VS_trial(trialNums, Y, rewardLabels, difficultyLabels, "Label", "EMG energy: "+muscleLabel(ch)+" (μV)", "OutputFolder", outputFolder+muscleLabel(ch), ...
            "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1700)
    % % EMG energy vs reward size
    % figure; hold on;
    % for j=1:ndifficulties
    %     Y = zeros(3,1);
    %     Yerr = zeros(3,1);
    %     for i=1:nrewards
    %         curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
    %         Y(i) = mean(energy(curInds, ch));
    %         Yerr(i) = std(energy(curInds, ch)) / sqrt(sum(curInds));
    %     end
    %     errorbar(1:nrewards, Y,Yerr, Color=diffColors(j, :), LineWidth=2);
    % end
    % set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    % xlim([0.7 3.3]); xticks(1:3); xticklabels(["S", "M", "L"]); ylabel("EMG energy: "+muscleLabel(ch)+" (μV)"); legend(["Tiny", "Huge"], Location="best");
    % set(gcf,'position',[0,0,550,550]);
    % saveas(gcf, outputFolder+muscleLabel(ch)+"-vs-Reward.jpg"); close all;

    % % EMG energy vs delay time
    % figure; hold on;

    % for j=1:ndifficulties
    %     for i=1:nrewards
    %         Y = zeros(nDelayTimes, 1); Yerr = zeros(nDelayTimes, 1);
    %         for k=1:nDelayTimes
    %             curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & delayTimes == DelayTimes(k);
    %             Y(k) = mean(energy(curInds, ch));
    %             Yerr(k) = std(energy(curInds, ch)) / sqrt(sum(curInds));
    %         end
    %         h(j) = errorbar(1:nDelayTimes, Y,Yerr, 'Color', rewColors(i, :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
    %     end
    % end
    % set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    % xlim([0.7 nDelayTimes+.3]); xticks(1:3:nDelayTimes); xticklabels(DelayTimes(1:3:end));
    % xlabel("Delay Length"); ylabel("EMG energy: "+muscleLabel(ch)+" (μV)"); legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
    % saveas(gcf, outputFolder+muscleLabel(ch)+"-vs-DelayLength.jpg"); close all;

    % % EMG energy vs direction
    % figure; hold on;
    % for j=1:ndifficulties
    %     for i=1:nrewards
    %         Y = zeros(ndirections, 1); Yerr = zeros(ndirections, 1);
    %         for k=1:ndirections
    %             curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & directionLabels == directions(k);
    %             Y(k) = mean(energy(curInds, ch));
    %             Yerr(k) = std(energy(curInds, ch)) / sqrt(sum(curInds));
    %         end
    %         h(j) = errorbar(1:ndirections, Y,Yerr, 'Color', rewColors(i, :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
    %     end
    % end
    % set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    % xlim([0.7 ndirections+.3]); xticks(1:ndirections); xticklabels([0 45 90 135 180 225 270 315]);
    % xlabel("Direction"); ylabel("EMG energy: "+muscleLabel(ch)+" (μV)"); legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
    % saveas(gcf, outputFolder+muscleLabel(ch)+"-vs-Direction.jpg"); close all;
end