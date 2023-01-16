function [trialData,taskInfo,neuralData] = pp_spikeAlignmentToTrial(taskAlignedSpiketimes,synchInfo,trialData,taskInfo)
% This function preprocesses the neural spike data by taking the full
% spiketrains of the data and properly apportioning them to the relevant
% trials.
%
% Inputs:
% - taskAlignedSpiketimes: structure with channel/sort labels, spiketimes,
% and mean waveform shapes for each unit
% - synchInfo: structure with many fields of synch information, but
% for this, all we use is the "taskSynchTrialTimes", which indicates when
% state 1 (center target appearance) first occurred for each trial ([ntrials
% x 1])
% - trialData: [ntrials x 1] structure with the trial parameter data and
% kinematics 
% - taskInfo: structure with various task parameters and information
%
% Outputs:
% - trialData: same as input, but with 1 more field:
% - neuralData: structure with the following fields:
%     - channel: [nunits x 1] channel # for each unit
%     - sort: [nunits x 1] sort # for each unit
%     - spikeMatrix: [nunits x ntime] boolean of if a spike occurred (at
%     ms precision)
% - taskInfo: same as input, but with 1 more field:
%   - meanWaveforms: [nunits x waveformLength] average spike waveform for
%   each unit
% Adam Smoulder, 11/2/22

disp('Aligning spike data to trial data')

% Get out unit info we'll need for every trial
unitLabels = taskAlignedSpiketimes.unitLabels;
meanWaveforms = taskAlignedSpiketimes.meanWaveforms;
spikeTimes = taskAlignedSpiketimes.spikeTimes;
nunits = length(unitLabels);
assert(length(meanWaveforms)==nunits)
assert(length(spikeTimes)==nunits)
channel = int32(cellfun(@(x) str2double(x(5:7)), unitLabels));
sortCodes = 'Uabcdef'; % unsorted will be 0, a will be 1, etc.
sort = int32(cellfun(@(x) find(sortCodes==x(end))-1, unitLabels));
waveformLength = median(cellfun(@(x) length(x), meanWaveforms)); % We assume at least half the units have waveforms
meanWaveformMatrix = nan(nunits,waveformLength); 
for j = 1:nunits
    if ~isempty(meanWaveforms{j})
        meanWaveformMatrix(j,:) = meanWaveforms{j};
    end
end; clear j

% We loop over each trial, get the relevant time window, get the spikes
% from each unit that are in that window, and proceed. 
ntrials = length(trialData);
synchTimes = synchInfo.taskSynchTrialTimes';
assert(length(synchTimes)==ntrials)
clear taskAlignedSpiketimes synchInfo
for i = 1:ntrials
    % Get the time for the trial and when the synch occurs (state 1)
    trialStartTime = double(trialData(i).time(1));
    trialEndTime = double(trialData(i).time(end));
    synchTimeInTrial = trialData(i).stateTable(2,find(trialData(i).stateTable(1,:)==1,1));
    time = trialData(i).time;
    ntime = length(time);
    
    % Find the spikes in this window
    spikeTimes_aligned = cellfun(@(x) floor((x-synchTimes(i))*1000)-synchTimeInTrial,spikeTimes,'uniformoutput',false); % align
    spikeTimes_aligned = cellfun(@(x) x(x >= trialStartTime & x <= trialEndTime),spikeTimes_aligned,'uniformoutput',false); % keep only ones in the trial
    
    % Make a spike matrix out of them
    spikeMatrix = false(nunits,ntime);
    for j = 1:nunits
        spikeMatrix(j,ismember(time,spikeTimes_aligned{j})) = true;
    end; clear j
    
    % Assign the relevant information
    neuralData(i).channel = channel;
    neuralData(i).sort = sort;
    neuralData(i).spikeMatrix = spikeMatrix;
    neuralData(i).meanWaveforms = meanWaveformMatrix;
%     trialData(i).neuralData = neuralData;
    clear spikeMatrix
    
    if mod(i,25)==0,disp(['Aligned spikes for trial ' num2str(i)]); end
end; clear i
taskInfo.ppInfo_NER.meanWaveforms = meanWaveformMatrix;

disp('Done aligning spikes to trials')

end

