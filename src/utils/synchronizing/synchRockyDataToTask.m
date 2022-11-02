%%% This script synchronizes Rocky's analog (.ns4) and neuron/threshold
%%% crossing data (.nex) to task data by identifying the timing of
%%% synchronization pulses.
%%%
%%% Synch pulses were sent from the task computer every time the center
%%% target appeared in a trial. These pulses were received by the Blackrock
%%% Cereplex box in analog input 16 as a downward pulse that lasted about
%%% 0.2ms. Each day of experiments from 02/2022 on should have at minimum 2
%%% minutes of no task at the start of each recording. From this, we
%%% identify the baseline voltage of analog input 16, then look for
%%% significant deviations from that to identify synch times. Because all
%%% other analog data was recorded through the same system at the same
%%% sampling frequency (10KHz), this effectively allows us to synchronize
%%% analog data to task events.
%%%
%%% Frustratingly, the system appears to have begun dropping synch pulses
%%% sometime between April and August 2022. Logic is in place to fill in
%%% the missing synch pulses using trial file times, though these should be
%%% visually validated (plots are in place for this purpose).
%%% 
%%% (Un-)Sorted neural data threshold crossing times (.nex) are effectively
%%% on the same clock, though importantly, the Blackrock system has an odd
%%% tendency to "reset" the time - typically at least once per session
%%% about 3-4 seconds in, then sporadically beyond that. We correct for
%%% the reset time in these different "chunks" to properly align spikes to
%%% the analog and task data.
%%%
%%% To perform synching, you'll need to provide the following:
%%% - name of anterior/posterior .nex files (if not using posterior, just
%%% pass in '')
%%% - name of anterior/posterior waveform file means (.txt; once again, if
%%% not posterior, just pass in '')
%%% - name of analog data file (.ns4)
%%% - folder names containing trial data for this task
%%%     - can take multiple folders for tasks split up over blocks that
%%%     were paused; pass each path in as an element of the cell array
%%% - the muscle names corresponding to each EMG channel (see expt log)
%%% - an output folder location
%%% - a folder path to fieldtrip's fileio folder (for loading .nex)
%%% - a folder path to NPMK (for loading .ns4)
%%%
%%% Tested successfully on datasets from the following:
%%% - original choking expt (3/3/22)
%%%     - No drops (or one?) - easy case
%%% - arousal interrupt choking expt, ANT only (3/25)
%%%     - No issues, occasionally 1 extra pulse
%%% - T/H data, ANT only (4/11/22)
%%%     - No issues, occasionally 1 extra pulse
%%% - Water/Juice data, ANT only (8/31/22)
%%%     - Tons of synch pulse drops; this seems to do well
%%%
%%% Tested unsuccessfully on:
%%% - Water/Juice day 8/30/22
%%%     - Data between trials 21-27 in block 3 have about 6.07s missing
%%%     from the analog data due to some weird system failure (see expt
%%%     log); synch pulses up to 21 and for 27 and beyond are present (27
%%%     and beyond needs the 6.07s offset), stuff inbetween is probs no
%%%     good
%%%
%%% Untested tasks: blocked reward range, RBD, focus
%%%
%%% Outputs data file with the following 3 structures:
%%% - analogData: structure with all analog data from .ns4 file. Contains the following fields:
%%%     - channelLabels: [nchannels x 1] cell array with the name for each of the analog data channels. Common channel label names include:
%%%         EMG_N (where N is a number, representing the Nth EMG channel), EKG, eyeX (gaze position X coordinate), eyeY (gaze position Y coordinate), 
%%%         eyeP (pupil size), HRM (not used - used to be for measuring heart rate with pulse oximeter), N/A (not used channel)
%%%     - EMGMuscleNames: [1 x nEMG] cell array with the names of the EMG muscles recorded corresponding to the EMG channels named EMG_1, EMG_2, etc.
%%%         Blank values indicate that no muscle was recorded for that channel
%%%     - fs: scalar, sampling frequency of the analog data (Hz). Should be 10000Hz for .ns4 files.
%%%     - time: [ntime x 1] double array with the aligned time values for the analog data (seconds).
%%%     - data: [ntime x nchannels] int16 array with the analog data. 
%%% - synchInfo: structure with information used from synchronization. Contains the following fields:
%%%     - taskSynchTrialTimes: [ntrials x 1] double array with the time of the center target appearance for each trial (seconds)
%%%     - synchUnitChan: string, name of analog input channel that had the synch pulse (this channel was removed after synchronization)
%%%     - sortedDataFile_ANT: string, name of sorted neural data file used for ANT array
%%%     - sortedDataFile_POST: string, name of sorted neural data file used for POST array
%%%     - BRFileFS: scalar, sampling rate of neural data recording (should be 30KHz)
%%%     - outputFolder: folder where I stored the synchronized file on my computer
%%% - taskAlignedSpiketimes: structure with neural data information. Contains the following fields:
%%%     - unitLabels: [nunits x 1] cell array with the name for each unit
%%%         Format: “Chan###x”, where ### = the channel number, and x = sort id (i.e. sort 1 of that channel is “a”, 2 is “b”) Channels 1-96 are the 
%%%         anterior array, channels 97-192 are the posterior array
%%%     - spikeTimes: [nunits x 1] cell array, where each cell contains a double array of the time of threshold crossing for that units’ sorted spikes
%%%     - meanWaveforms: [nunit x 1] cell array, where each cell has the mean waveform shape for that unit
%%%
%%% I recommend running this section by section for each day to make sure
%%% eveyrthing goes smoothly.
%%%
%%% Adam Smoulder, last edited 11/2/22

  
% ATTEMPT FOR OLD CHOKING; no issues
sortedDataFilename_ANT = 'F:\~~~RockyData\BORTRig\Blackrock_anterior\20220216to20220303_chokingReachingExperimentsAndChoice\Rocky_Anterior_20220302_delayedCOut_sorted.nex';
waveformMeanFilename_ANT = 'F:\~~~RockyData\BORTRig\Blackrock_anterior\20220216to20220303_chokingReachingExperimentsAndChoice\Rocky_Anterior_20220302_delayedCOut_sorted.txt';
sortedDataFilename_POST = 'F:\~~~RockyData\BORTRig\Blackrock_posterior\20220216to20220303_chokingReachingExperimentsAndChoice\Rocky_Posterior_20220302_delayedCOut_sorted.nex';
waveformMeanFilename_POST = 'F:\~~~RockyData\BORTRig\Blackrock_posterior\20220216to20220303_chokingReachingExperimentsAndChoice\Rocky_Posterior_20220302_delayedCOut_sorted.txt';
analogDataFilename = 'F:\~~~RockyData\BORTRig\Blackrock_anterior\20220216to20220303_chokingReachingExperimentsAndChoice\Rocky_Anterior_20220302_delayedCOut.ns4';
taskDataFolders = {... % don't forget the \ at the end
    'F:\~~~RockyData\BORTRig\TaskData\20220216to20220303_chokingReachingExperimentsAndChoice\centerOutData_20220302_153642\'
    };
outputFolder = '..\..\..\data\synchronized\'; % where to save the file; by default, it will use the datestring from the task.
EMGMuscleNames = {'Tric','PDel','Trap','ADel','LBic','',''};





disp('Synchronizing and aligning data with no SPM')

% Parameters that should not change across sessions
analogDataLabelMapping = {... % Analog channel names
    ['ainp1' char(0)],'EMG_1';... % for some god-forsaken reason, MATLAB's space is a different character from theirs??? ugh...
    ['ainp2' char(0)],'EMG_2';...
    ['ainp3' char(0)],'EMG_3';...
    ['ainp4' char(0)],'EMG_4';...
    ['ainp5' char(0)],'EMG_5';...
    ['ainp6' char(0)],'EMG_6';...
    ['ainp7' char(0)],'EMG_7';...
    ['ainp8' char(0)],'EKG';...
    ['ainp9' char(0)],'eyeX';...
    'ainp10','eyeY';...
    'ainp11','eyeP';...
    'ainp12','HRM';...
    'ainp13','N/A';...
    'ainp14','N/A';...
    'ainp15','N/A';...
    'ainp16','Synch';...
    };

% Parameters from recording and file naming conventions; shouldn't change
BRFileFS = 30000;                   % BR machine samples at 30KHz
fieldtripPath = 'D:\fieldtrip-20160712\fieldtrip-20160712\fileio\';  % we need fileio from fieldtrip to load the .nex files
npmkPath = 'F:\~~~RockyData\BORTRig\~NPMK'; % we need this for loading the analog data file
synchChanName = 'ainp16';      % which analog signal is the synch pulse?

assert(~strcmp(sortedDataFilename_ANT(end-3:end),'.plx'),'Do not use PLX! Takes too much memory (loads all waveforms). Convert to NEX')


%% Add paths and load analog data
disp('Adding paths and loading data')
addpath(fieldtripPath)
addpath(genpath(npmkPath))
addpath(pwd)

% Load analog data
analogStruct = openNSx(analogDataFilename);
disp('Data loaded')


%% Extract info we need
Fs = analogStruct.MetaTags.SamplingFreq;
nchans = analogStruct.MetaTags.ChannelCount;
chanLabels = cell(nchans,1);
for c = 1:nchans
    chanLabels{c} = analogStruct.ElectrodesInfo(c).Label;
    curChanData = cellArrayToVector(cellfun(@(x) x(c,:), analogStruct.Data,'uniformoutput',false)); % UNITS ARE 1/4s OF uVs! so this /4 = uV
    if c == 1
        analogDataTraces = curChanData;
    else
        analogDataTraces(:,end+1) = curChanData;
    end
    pause(1) % no idea why, but this stops my memory from blowing up and cpu from crashing, so pausing it is
    disp(['Extracted ch' num2str(c)]) 
end; clear c
dt = 1/Fs;
ntime_NSx = size(analogDataTraces,1);
time_NSx = (0:dt:(ntime_NSx-1)*dt)'; % actual value isn't as big of a deal as step size being right

% Map the channel names to their meaning
newLabels = cell(nchans,1);
for c = 1:nchans
    newLabelInd = find(cellfun(@(x) contains(chanLabels{c},x), analogDataLabelMapping(:,1)),1,'first');
    newLabels{c} = analogDataLabelMapping{newLabelInd,2};
end; clear c
analogData.channelLabels = newLabels;
analogData.EMGMuscleNames = EMGMuscleNames;
analogData.fs = Fs;
analogData.time = time_NSx;
analogData.data = analogDataTraces; % this is the original data precision and stores 1/4 the size of double


%% Threshold synch signals
disp('Extracting synch signals; assuming first 30s have no synch pulses')
synchChanInd = find(cellfun(@(x) contains(x,synchChanName), chanLabels));
if ~isempty(synchChanInd)
    synchTrace_NSx = double(analogDataTraces(:,synchChanInd));
    analogDataTraces(:,synchChanInd) = [];
    chanLabels(synchChanInd) = [];
    nchans = nchans-1;
end
synchInds_NSx  = synchTrace_NSx < median(synchTrace_NSx)-20*std(synchTrace_NSx(1:30*Fs)); 

% Synch is above threshold multiple indices in a row, average them
synchTimes_NSx = time_NSx(synchInds_NSx);
multipleTrigInds = find(diff(synchTimes_NSx)-dt < 0.1); % no two synchs should be < 0.1s apart
synchTimes_NSx_old = synchTimes_NSx;
for i = 1:length(multipleTrigInds)
    curSynchTime = synchTimes_NSx(multipleTrigInds(i));
    curMultTrigInds = find(abs(curSynchTime-synchTimes_NSx)-dt < 0.1);
    synchTimes_NSx(curMultTrigInds) = mean(synchTimes_NSx_old(curMultTrigInds));
end; clear i
synchTimes_NSx(multipleTrigInds) = [];
multipleTrigInds = find([1; diff(synchTimes_NSx)-dt] < 1E-10);
assert(isempty(multipleTrigInds))


%% Extract info from task data for synchronization
disp('Extracting synch info from trial data')
ntrials = 0;
trialDurations = nan(5000,1);
trialSynchTimes = nan(5000,1);
fileCreationTimes = nan(5000,1);
dates = cell(5000,1); % used for the date string
counter = 1;
for f = 1:length(taskDataFolders)
    curFiles = dir(taskDataFolders{f});
    for i = 1:length(curFiles)
        curFname = [taskDataFolders{f} curFiles(i).name];
        if length(curFname)-length(taskDataFolders{f}) > 6 && strcmp(curFname((1:6)+length(taskDataFolders{f})),'Trial0')
            ntrials = ntrials+1;
            curTrial = load(curFname);
            trialDurations(counter) = double(curTrial.trialData.stateTable(end,end))/1000; % last value is last timept in ms
            trialSynchTimes(counter) = double(curTrial.trialData.stateTable(2,1))/1000; % first value is synch time (center onset)
            fileCreationTimes(counter) = curFiles(i).datenum; % Sadly precision is like 2s, so we'll have to get creative if needed
            dates{counter} = curFiles(i).date;
            counter = counter+1;
        end
    end; clear i
end; clear f
trialDurations(counter:end) = []; trialSynchTimes(counter:end) = []; fileCreationTimes(counter:end) = []; dates(counter:end) = [];
fileCreationTimes = fileCreationTimes*24*60*60; % convert to s


%% Now to figure out where each trial lands in time_NSx...

% Get estimates of synch time in cpu framing from file creation times
fcts = fileCreationTimes-trialDurations+trialSynchTimes;
fcts = fcts-min(fcts);

% Evaluate the current case...
% ...if numbers match and correlation is high, seems like we're clear
if length(fcts)==length(synchTimes_NSx) && max(abs(diff(synchTimes_NSx)-diff(fcts))) < 2.1 % should be < 2s error as that's the file creation precision
    disp('All synchs accounted for')
    synchTimes_NSx_final = synchTimes_NSx; % everything is going to be in analog data's time
    %     synchInfo.badSynchInds = badInds; % these trials cannot be properly synchronized and should be excluded
else %... but if not, we have to figure out which pulses map to which trials
    if length(synchTimes_NSx) > length(fcts) % easier case; we have more synch pulses than trials = probs just one extra or something
        if length(synchTimes_NSx)-length(fcts)==1 && max(abs(diff(synchTimes_NSx(1:end-1))-diff(fcts))) < 2.1 % just 1 extra pulse tacked-on; easy case
            disp('All synchs accounted for (one extra pulse for last (unsaved) trial)')
            synchTimes_NSx_final = synchTimes_NSx(1:end-1); % everything is going to be in analog data's time
    %     synchInfo.badSynchInds = badInds; % these trials cannot be properly synchronized and should be excluded
        else
            error('#pulses > #trials by more than 1, investigate')
        end
    else % more trials than pulses = we have drops...
        disp('Synch mismatch; finding best mapping of synch pulses to trials')
        synchTimes_NSx_working = synchTimes_NSx; % what we'll use as we work
        
        % First, identify when long pauses occur in both NSx and task files
        longPauseMinLength = 30; % 30 seems to work most days?
        diffs_fct = diff(fcts);
        longPauses_fct = find([0;diffs_fct>longPauseMinLength]); % longer than 30s b/w file creation = long pause/block switch
        diffs_nsx = diff(synchTimes_NSx_working);
        longPauses_nsx = find([0;diffs_nsx>longPauseMinLength]);
        if length(longPauses_nsx) == length(longPauses_fct)+1 % probably warmup block; remove these pulses
            disp('Removing warmup block pulses')
            synchTimes_NSx_working(1:longPauses_nsx(1)) = [];
            diffs_nsx = diff(synchTimes_NSx_working);
            longPauses_nsx = find([0;diffs_nsx>longPauseMinLength]);
        end
        assert(length(longPauses_nsx)==length(longPauses_fct),'Different number of long pauses; investigate')
        
        % Split data up by the long pauses and figure out stuff
        synchTimes_NSx_final = nan(size(fcts));
        errors_final = nan(size(fcts));
        for b = 1:length(longPauses_fct)+1
            disp(['Aligning synchs for block ' num2str(b)])
            if b == 1 % get trial indices for the current pause
                inds_fct = 1:longPauses_fct(b)-1;
                inds_nsx = 1:longPauses_nsx(b)-1;
            elseif b == length(longPauses_fct)+1
                inds_fct =  longPauses_fct(b-1):length(fcts);
                inds_nsx = longPauses_nsx(b-1):length(synchTimes_NSx_working);
            else
                inds_fct = longPauses_fct(b-1):longPauses_fct(b)-1;
                inds_nsx = longPauses_nsx(b-1):longPauses_nsx(b)-1;
            end
            nsx_b = synchTimes_NSx_working(inds_nsx);
            nfct_b = length(inds_fct); % # trials in this block
            nnsx_b = length(nsx_b); % # synch pulses in this block
            
            % Trial durations (shifted by trial synch) should match up with good precision (<< 2s) to synch pulses; use this to find mapping
            % As a result, pairwise pulse vs. synch time error should show obvious structure to indicate which trial each pulse corrseponds to
            diff_synchTimes_b = trialDurations(inds_fct(1:end-1))-trialSynchTimes(inds_fct(1:end-1))+trialSynchTimes(inds_fct(2:end))-0.001; % ~1ms is the avg lag
            sts_b = [0 ; cumsum(diff_synchTimes_b)]; % synch times
            refPt = 1; % assumes 1st pulse references trial #<refPt> (starts with 1; if error is high, we step up, assuming early pulses are dropped)
            while refPt==1 || mean(abs(errors)) > 0.010
                if refPt~=1
                    disp(['High error in synch time reconstruction; trying assuming ' num2str(refPt-1) ' trial pulse-drop'])
                end
                errorMat = (sts_b-sts_b(refPt))-(nsx_b'-nsx_b(1));
                [errors,inds] = min(abs(errorMat));
                for i = 1:nnsx_b
                    errors(i) = errorMat(inds(i),i);
                end; clear i
                figure; subplot(1,2,1);
                imagesc(abs(errorMat)); colorbar; set(gca,'clim',[0 0.05]); title([num2str(refPt-1) ' synch pulse(s) assumed dropped, Block ' num2str(b) ', offset = ' num2str(0)])
                subplot(1,2,2); plot(errors,'o'); title(['Avg abs error = ' num2str(mean(abs(errors)))])
                set(gcf,'position',[-3839 -222 1920 963])
                pause(1)
                offset = median(errors);
                errorMat = (sts_b-sts_b(refPt)-offset)-(nsx_b'-nsx_b(1));
                [errors,inds] = min(abs(errorMat));
                for i = 1:nnsx_b
                    errors(i) = errorMat(inds(i),i);
                end; clear i
                figure; subplot(1,2,1);
                imagesc(abs(errorMat)); colorbar; set(gca,'clim',[0 0.05]); title([num2str(refPt-1) ' synch pulse(s) assumed dropped, Block ' num2str(b) ', offset = ' num2str(offset)])
                subplot(1,2,2); plot(errors,'o'); title(['Avg abs error = ' num2str(mean(abs(errors)))])
                set(gcf,'position',[-3839 -222 1920 963])
                pause(1)
                
                refPt = refPt+1;
            end
            
            % Assign the appropriate trials to the appropriate indices
            nsx_b_final = nan(nfct_b,1);
            for i = 1:nnsx_b
                if abs(errors(i)) < 0.010 % less than 10ms is pretty good
                    nsx_b_final(inds(i)) = nsx_b(i);
                end
            end; clear i
            missingPulseInds = find(isnan(nsx_b_final));
            
            % Fill in the missing pulse times. We do this by:
            % - get diff in estimated synch times from previous to missing pulse
            % - add this diff to the previous pulse time
            % If the first few pulses are dropped, we have to work backward for those
            if missingPulseInds(1)==1 % work backward for first pulses dropped if needed
                firstNonmissing = find(~ismember(1:nfct_b,missingPulseInds),1);
                for p = firstNonmissing-1:-1:1
                    curRevDiff = sts_b(p+1)-sts_b(p);
                    nsx_b_final(p) = nsx_b_final(p+1)-curRevDiff;
                end; clear p
            end
            missingPulseInds = find(isnan(nsx_b_final));
            for p = 1:length(missingPulseInds) % o/w work forward
                curDiff = diff(sts_b(missingPulseInds(p)+[-1 0]));
                nsx_b_final(missingPulseInds(p)) = nsx_b_final(missingPulseInds(p)-1)+curDiff;
            end; clear p
            errors_b = (nsx_b_final-nsx_b_final(1))-(sts_b-sts_b(1));
            offset = mean(errors_b);
            errors_b = errors_b - offset;
            
            % Plot the errors and move onto the next block
            figure
            plot(errors_b,'o'); title('Est. synch time to pulse errors; look for outliers')
            set(gcf,'position',[-3158 154 560 420])
            synchTimes_NSx_final(inds_fct) = nsx_b_final;
            errors_final(inds_fct) = errors_b;
            pause(1)
        end; clear b
    end
end

synchInfo.taskSynchTrialTimes = synchTimes_NSx_final; % everything is going to be in analog data's time
% synchInfo.badSynchInds = badInds; % these trials cannot be properly synchronized and should be excluded for analyses
synchInfo.dates = dates;
synchInfo.synchUnitChan = synchChanName;
synchInfo.sortedDataFile_ANT = sortedDataFilename_ANT;
synchInfo.sortedDataFile_POST = sortedDataFilename_POST;
synchInfo.BRFileFS = BRFileFS;
synchInfo.outputFolder = outputFolder;


%% Clear some space!
clear synchTrace_NSx synchInds_NSx curChanData analogDataTraces 
analogStruct = rmfield(analogStruct,'Data');


%% Load offline sorted data
ANT_sortedSpikesStructure = ft_read_spike(sortedDataFilename_ANT);
ANT_spikeTimes_sorted = cellfun(@double,ANT_sortedSpikesStructure.timestamp,'uni',0)';
ANT_spikeTimes_sorted = cellfun(@(x) x/BRFileFS, ANT_spikeTimes_sorted,'uniformoutput',false); % convert from samples to s
disp('Loaded anterior sorted spikes')

if ~isempty(sortedDataFilename_POST)
    POST_sortedSpikesStructure = ft_read_spike(sortedDataFilename_POST);
    POST_spikeTimes_sorted = cellfun(@double,POST_sortedSpikesStructure.timestamp,'uni',0)';
    POST_spikeTimes_sorted = cellfun(@(x) x/BRFileFS, POST_spikeTimes_sorted,'uniformoutput',false); % convert from samples to s
    disp('Loaded posterior sorted spikes')
else
    disp('No posterior array data used')
end

% Get waveform means
sortCodes = ['Uabcde']; % u = 0, a = 1, b = 2, etc.
if ~isempty(waveformMeanFilename_ANT)
    ANT_waveformMeans = cell(length(ANT_sortedSpikesStructure.label),1);
    txtFileData = table2array(readtable(waveformMeanFilename_ANT));
    channel = txtFileData(:,1);
    sort = txtFileData(:,2);
    meanWaveformData = txtFileData(:,3:end);
    for i = 1:length(ANT_sortedSpikesStructure.label) % for each unit
        curUnitLabel = ANT_sortedSpikesStructure.label{i}; % "Chan###L" is the format; ### = channel #, L = sort code label
        curChan = str2double(curUnitLabel(end-3:end-1));
        curUnitNum = find(sortCodes==curUnitLabel(end))-1; % 0 = U, 1 = a, etc.
        if ~isempty(curUnitNum) % if unit is found in mean waveform file
            curUnitIndex = find(ismember([channel sort],[curChan curUnitNum],'rows'),1);
            ANT_waveformMeans{i} = meanWaveformData(curUnitIndex,:)';
        end
    end; clear i
end
if ~isempty(waveformMeanFilename_POST)
    POST_waveformMeans = cell(length(POST_sortedSpikesStructure.label),1);
    txtFileData = table2array(readtable(waveformMeanFilename_POST));
    channel = txtFileData(:,1);
    sort = txtFileData(:,2);
    meanWaveformData = txtFileData(:,3:end);
    for i = 1:length(POST_sortedSpikesStructure.label) % for each unit
        curUnitLabel = POST_sortedSpikesStructure.label{i}; % "Chan###L" is the format; ### = channel #, L = sort code label
        curChan = str2double(curUnitLabel(end-3:end-1));
        curUnitNum = find(sortCodes==curUnitLabel(end))-1; % 0 = U, 1 = a, etc.
        curUnitIndex = find(ismember([channel sort],[curChan curUnitNum],'rows'),1);
        if ~isempty(curUnitIndex) % if unit is found in mean waveform file
            POST_waveformMeans{i} = meanWaveformData(curUnitIndex,:)';
        end
    end; clear i
end
disp('Loaded waveform mean files')

% We need to update all of the poterior array channel names before
% concatentation (add 96). Format is Chan###L, where ### is the channel
% number (always 3 digits) and L is the unit label (u,a,b,c,d,e)
if ~isempty(sortedDataFilename_POST)
    for i = 1:length(POST_sortedSpikesStructure.label)
        oldLabel = POST_sortedSpikesStructure.label{i};
        POST_sortedSpikesStructure.label{i} = [...
            oldLabel(1:4),... % 'Chan'
            sprintf('%03d',str2double(oldLabel(5:7))+96),... % channel ### (96 added)
            oldLabel(end),... % the unit label
            ];
    end; clear i
end

% Identify unsorted units to remove
if strcmp(waveformMeanFilename_ANT,'fake96chanUnsortedWaveformMeans.txt')
    disp('USING UNSORTED DATA')
    ANT_unsortedUnits = false(96,1);
    POST_unsortedUnits = false(96,1);
else
    disp('Labeling unsorted units to be removed')
    ANT_unsortedUnits = cellfun(@(x) any(ismember(x,'U')),ANT_sortedSpikesStructure.label);
    if ~isempty(sortedDataFilename_POST)
        POST_unsortedUnits = cellfun(@(x) any(ismember(x,'U')),POST_sortedSpikesStructure.label);
    end
end

% Combine labels
if ~isempty(sortedDataFilename_POST)
    taskAlignedSpiketimes.unitLabels = [...
        ANT_sortedSpikesStructure.label(~ANT_unsortedUnits)';
        POST_sortedSpikesStructure.label(~POST_unsortedUnits)'...
        ];
    taskAlignedSpiketimes.meanWaveforms = [...
        ANT_waveformMeans(~ANT_unsortedUnits) ; ...
        POST_waveformMeans(~POST_unsortedUnits)
        ];
    spikeTimes_sorted = [ANT_spikeTimes_sorted(~ANT_unsortedUnits) ; POST_spikeTimes_sorted(~POST_unsortedUnits)];
else
    taskAlignedSpiketimes.unitLabels = ANT_sortedSpikesStructure.label(~ANT_unsortedUnits)';
    taskAlignedSpiketimes.meanWaveforms = ANT_waveformMeans(~ANT_unsortedUnits);
    spikeTimes_sorted = ANT_spikeTimes_sorted(~ANT_unsortedUnits);
end


%% Align spikes with analog data. Every time the Blackrock system chunks 
%  data, the clock resets to 0, so we'll have to adjust the spiketimes
%  accordingly. Mindbogglingly frustratingly, only in certain datasets does
%  a chunk split reset the clock! Check these results closely each time...
%
%  See comments at the bottom of the script for details.

% Identify if clock resets occurred for each of the chunks
chunkDurations = analogStruct.MetaTags.DataPointsSec;
chunkCumDurations = cumsum(chunkDurations);
nchunks = length(chunkDurations);
spikeTimes_sorted_orig = spikeTimes_sorted;
nunits = length(spikeTimes_sorted_orig);
clockResetInds = cellfun(@(x) find(diff(x)<0)+1, spikeTimes_sorted,'uniformoutput',false);
nresets = max(cellfun(@(x) length(x), clockResetInds));
if nresets==0 % do nothing
    chunkClockReset = false(nchunks-1,1);
elseif nresets==nchunks-1 % time drop for each chunk
    chunkClockReset = true(nchunks-1,1);
else % time drops, but not for each chunk...need to identify which
    validUnits = cellfun(@(x) length(x), clockResetInds)==nresets; % only check with units that show all drops
    chunkClockReset = false(nchunks-1,1);
    chunkCumDurations_working = chunkCumDurations;
    for r = 1:nresets
        % Estimate time of drop based on latest spike
        if r == 1
            maxPreDropSpiketimes = cellfun(@(x,inds) max(x(1:inds(r))),spikeTimes_sorted(validUnits),clockResetInds(validUnits));
        else
            maxPreDropSpiketimes = cellfun(@(x,inds) max(x(inds(r-1):inds(r))),spikeTimes_sorted(validUnits),clockResetInds(validUnits));
        end
        estResetTime = max(maxPreDropSpiketimes); % should correspond with a chunk or sum of a few cumulative ones
        
        % Find which cumulative sum of chunks best matches this length
        [error,chunkInd] = min(abs(chunkCumDurations_working-estResetTime));
        chunkClockReset(chunkInd) = true;
        
        % Update cumulative chunk lengths to match clock reset
        chunkCumDurations_working = chunkCumDurations_working-chunkCumDurations_working(chunkInd);
    end; clear r
end
chunkClockResetInds = find(chunkClockReset); % make it indices instead
clockResetTimes = chunkCumDurations(chunkClockResetInds);
updateValues = clockResetTimes-[0 clockResetTimes(1:end-1)]; % value to add to compensate for drop

% Now update spikes for when clock resets occurred
for u = 1:nunits
    if isempty(clockResetInds{u})% no resets = assume first spike happens after any drops; see comments at bottom of code for details
        spikeTimes_sorted{u} = spikeTimes_sorted{u}+sum(updateValues);
    else
        % We identify which reset each of the ones for this unit corresponds to,
        % then add all preceding drops that haven't yet been accounted for.
        updateValTotal = 0;
        for r = 1:length(clockResetInds{u})
            estResetTime = spikeTimes_sorted{u}(clockResetInds{u}(r)-1)-updateValTotal;
            [error,curResetInd] = min(abs(estResetTime-(clockResetTimes-updateValTotal)));
            curUpdateVal = clockResetTimes(curResetInd)-updateValTotal;
            spikeTimes_sorted{u}(clockResetInds{u}(r):end) = ...
                spikeTimes_sorted{u}(clockResetInds{u}(r):end)+curUpdateVal;
            updateValTotal = updateValTotal+curUpdateVal;
        end; clear r
    end
end; clear u
assert(~any(cellfun(@(x) any(diff(x)<0), spikeTimes_sorted)),'Fixing BR chunks did not work; investigate')

% Plot an example channel to check
[~,exampleInd] = max(cellfun(@(x) length(x), spikeTimes_sorted));
figure; hold on
plot(spikeTimes_sorted_orig{exampleInd},'o')
plot(spikeTimes_sorted{exampleInd},'-','linewidth',1)
legend('orig','updated')
title(['Neuron ' num2str(exampleInd) ' spiketimes; zoom in early and around breaks'])
set(gcf,'position',[-3839 -222 1920 963])
pause(1)
% To be clear, what this show is the original spike times for this neuron
% (blue circles) and the ones corrected for clock resets (red line). The
% red line should look nearly continuous and that it corrects for resets in
% the blue.

% Assign it
taskAlignedSpiketimes.spikeTimes = spikeTimes_sorted;


%% Finally, save!!!

% Get date number
dateStartInd = strfind(sortedDataFilename_ANT,'202'); % for 2021/2022
dateStartInd = dateStartInd(end); % in case there are multiple 202's in the path name...
dateString = sortedDataFilename_ANT(dateStartInd:dateStartInd+7)

outputFilename = [outputFolder 'Rocky_synchedSpikeAndAnalogData_noSPM_' dateString]
save(outputFilename,'taskAlignedSpiketimes','synchInfo','analogData','-v7.3')


disp('Done Saving!')


%% 
%%% FOR UPDATING SPIKES WHEN CLOCK RESETS OCCUR
%%% The basic gist is as follows:
%%% - BR data seems to be saved in data "chunks" of variable length
%%%     - Seemingly EVERY session has at least 2, where the chunks split
%%%     about 3-4 seconds in when the two recording boxes synchronize
%%% - Sometimes (but not always), when a new chunk begins, it's accompanied
%%% with a clock reset - that is, the time of all threshold crossings is
%%% reset to 0
%%%     - BUT conveniently, the indices of data sorted remain in order;
%%%     hence if a neuron had spikes at times [2 3 7 8] and a clock reset
%%%     occurred at 6s, it would look like [2 3 1 2]
%%% - So, assuming at least one neuron spikes fast enough to show any given
%%% clock reset, we find which chunk splits correspond with clock resets by
%%% seeing if any neuron has a diff in spike times < 0
%%% - We then update the spikes that are after the clock reset to add the
%%% appropriate chunk length
%%% - There are a few (rare) cases where it is ambiguous if a spike occurs
%%% before or after a clock reset.
%%%     - i.e., raw spikes at [2 6 8...] and we know a clock reset occurred
%%%     at 3. Is the 2 from before, or after the clock reset?
%%%     - These cases only occur in the case of (1) Very short chunks, and
%%%     (2) Very sparsely spiking neurons. The only time this would be
%%%     relevant to the experiment is if many VERY short blocks (i.e. < 1s
%%%     long) happened concecutively during the experiment
%%%     - Because of this, we just assume all spikes happened AFTER the
%%%     clock reset and correct time accordingly










