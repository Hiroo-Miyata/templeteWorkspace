%%% This script takes synchronized data and preprocesses all data streams.
%%% It saves 8 separate files for a session:
%%% - TSK: task parameter data (e.g., dir/rew labels, state trans., etc.)
%%% - KIN: kinematic data for hand position
%%% - EMG: muscle data
%%% - ECG: heart rate information
%%% - EYE: eye gaze and pupil data
%%% - NER: neural data
%%% - BEH: single trial behavioral metrics (e.g., peak speed, saccade RT)
%%% - taskInfo: metrics relevant across the entire session (i.e., task
%%% parameters, mean waveforms, filter information)
%%%
%%% Aside from taskInfo, all of these are structures with length ntrials
%%% that have fields with values corresponding to each trial. Aside from
%%% taskInfo/BEH, all field values within a trial have length ntime, where
%%% ntime is the number of milliseconds in that trial.
%%%
%%% It's recommended you run each of the preprocessing steps in sections to
%%% evaluate output as you go.
%%%
%%% Adam Smoulder, 11/2/2022

dbstop if error

% Edit folder/filenames and file tags (what they start with) here
spikeAndAnalogDataFile = ...
    'F:\~~~RockyData\BORTRig\~alignedSpikeAndAnalogData_synchDrops\Rocky_synchedSpikeAndAnalogData_noSPM_20220302.mat';
taskDataFolderList = {...
    'F:\~~~RockyData\BORTRig\TaskData\20220216to20220303_chokingReachingExperimentsAndChoice\centerOutData_20220302_153642\';
    };
functionFolder = ['..\']; % This is where the functions and classes needed are, at least beneath it
outputFolder = ['..\..\..\data\preprocessed\'] ; % where to save files



% Things that shouldn't change over days
trialDataTag = 'Trial0'; % what each trial data file starts with
rewardNames = {'Small','Medium','Large','Jackpot'}; % Searches target names for these

% Go to working folder
addpath(functionFolder)

% % Adam uses this to remove paths to other stuff
% rmpath(genpath('D:\AdamMatlab'))
% addpath(genpath(pwd))

%% Load the data
disp('Loading data')
dateString = taskDataFolderList{1}(end-15:end-1) % show it for confirmation
load([taskDataFolderList{1} 'taskInfo.mat'])
load(spikeAndAnalogDataFile)
disp('Loaded synchronized data')

%% Do an initial processing of the behavior and set up the file name to save
[trialData,taskInfo,kinematicData] = pp_preprocessTaskAndKinematicData(taskDataFolderList,rewardNames,synchInfo);
procTimeString = datestr(now,'yyyymmdd_HHMMSS');
saveFilenameBase = ['Rocky' ...
    taskDataFolderList{1}(end-15:end-8) ...
    '_Trials' num2str(trialData(1).trial) ...
    '_' num2str(trialData(end).trial) '_preprocessed_' ...
    ];
save([outputFolder saveFilenameBase 'TSK_' procTimeString],'trialData')
save([outputFolder saveFilenameBase 'KIN_' procTimeString],'kinematicData')
disp('Saved task and kinematic data')

%% Attach spikes to trialData
[trialData,taskInfo,neuralData] = pp_spikeAlignmentToTrial(taskAlignedSpiketimes,synchInfo,trialData,taskInfo);
save([outputFolder saveFilenameBase 'NER_' procTimeString],'neuralData')
disp('Saved neural data')

%% Preprocess EMG and ECG data
% TO-DO: put in EMG processing here
outputFormat = ["jpg", "svg", "fig"]; %"jpg", "svg", "fig"
plotoutputFolder = [outputFolder 'emgConfirmPlot/'];
[EMGData, EMGMetrics, ECGData] = emgPreprocess(trialData, analogData, outputFormat, plotoutputFolder);
save([outputFolder saveFilenameBase 'EMG_' procTimeString], 'EMGData', 'EMGMetrics')
save([outputFolder saveFilenameBase 'ECG_' procTimeString], 'ECGData')

%% Preprocess eye data
% TO-DO: put in Eye data processing here
% save([outputFolder saveFilenameBase 'EYE_' procTimeString],'eyeData')

%% Assess various data streams for behavioral metrics (i.e., peak speed,
% saccade reaction time, average heart rate, etc.)
% TO-DO: put in behavioral processing here
% save([outputFolder saveFilenameBase 'BEH_' procTimeString],'behaviorData')

disp('Completed preprocessing')





